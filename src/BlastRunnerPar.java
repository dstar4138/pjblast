/**
 * Parallel Version of BLAST.
 * Runs a search over a database list in parallel, searching each one on a separate process and reporting 
 * back the alignment ranges are. Also the parallel BLAST function makes use of a multicore environment to 
 * handle long sub word searches within the database entry itself.
 **/

import edu.rit.compbio.seq.Alignment;
import edu.rit.compbio.seq.AlignmentPrinter;
import edu.rit.compbio.seq.DefaultAlignmentStats;
import edu.rit.compbio.seq.ProteinDatabase;
import edu.rit.compbio.seq.ProteinSequence;
import edu.rit.compbio.seq.Sequence;

import edu.rit.pj.Comm;
import edu.rit.pj.WorkerTeam;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerLongForLoop;
import edu.rit.mp.buf.ObjectItemBuf;

import java.io.File;
import java.util.Hashtable;

public class BlastRunnerPar {

    public static void main( String[] args ) throws Exception{
        Comm.init(args);
        final Comm world = Comm.world();
        final int rank = world.rank();
        final int size = world.size();
         
        if( args.length != 5 ){ usage(); }
        // TODO: ProteinDatabase wont necessarily work if the user passes in a nucleotide file in.
        final ProteinDatabase pd = new ProteinDatabase( new File( args[1] ) , new File(args[2]) );
        final float percentage = (new Float(args[3])).floatValue();
        
        final Sequence query;
        final BLAST aligner;
        if(args[0].toLowerCase().startsWith("p")){
            query = new ProteinSequence(">Search Query", args[4]);
            aligner = new BLASTP();
        }else{
            query = new NucleotideSequence(">Search Query", args[4]);
            aligner = new BLASTN();
        }
        
        long t1 = System.currentTimeMillis();
        
        // Loop over database and find alignments.
        final Hashtable<Sequence, Alignment[]> processData = new Hashtable<Sequence, Alignment[]>();
        new WorkerTeam().execute( new WorkerRegion(){
            public void run() throws Exception {
                execute( 0L,  (long) (pd.getProteinCount()*percentage), new WorkerLongForLoop(){
                    public void run( long first, long last ){
                    	AlignmentPrinter out;
                        for( long i=first; i<=last; ++i){
                        	try{
	                            Alignment[] tmp = aligner.align( query, pd.getProteinSequence( i ) );
	                            
	                            // If we are process 0, print as we find alignments. If we are any other process,
	                            // hold our values until we finish all of them and then send it as a lump sum to
	                            // process 0. (This conserves our message calls and reduces latency).
	                            if(rank == 0){
	                            	out = new AlignmentPrinter(System.out,new DefaultAlignmentStats(pd.getProteinSequence(i).length()));
	                    			for(Alignment align : tmp)
	                       	 			out.printDetails(align, query, pd.getProteinSequence( i ));
	                            }else{
	                            	processData.put(pd.getProteinSequence( i ), tmp);
	                            }
                        	}catch( Exception e ){
                        		//LATER: Logging?
                        	}
                        }
                    }
                    
                    public void finish() throws Exception{
                    	ObjectItemBuf<Hashtable<Sequence, Alignment[]>> buf = new ObjectItemBuf<Hashtable<Sequence, Alignment[]>>();
                    	buf.item = processData;
                    	world.send(0, buf);
                    }
               });
             }
        });
        
        // Listen for alignments from other processes when they are done.
        if( rank == 0 ){
        	AlignmentPrinter out;
        	for(int i=1; i<=size; i++){
        		ObjectItemBuf<Hashtable<Sequence, Alignment[]>> buf = new ObjectItemBuf<Hashtable<Sequence, Alignment[]>>();
        		world.receive(i, buf); //TODO: validate CommStatus return?
        		
        		//Print sequences that this thread went over.
        		for(Sequence subject : buf.item.keySet()){
        			out = new AlignmentPrinter(System.out,new DefaultAlignmentStats(subject.length()));
        			for(Alignment align : buf.item.get(subject))
           	 			out.printDetails(align, query, subject);
        		}
        	}
        }
        
        long t2 = System.currentTimeMillis();
        System.out.printf("Running Time: %d msec%n", t2-t1);
        
    }

    private static void usage(){
        System.out.println( "Usage: java BlastRunnerSeq [p/n] <input database> <index file> <search coverage percentage> <query seqence>\n"+
                            "\t[p/n] - p if proteins, n if nucleotides.\n"+
                            "\t<input database> - Path to a file with FASTA syle sequences.\n" +
                            "\t<database indexs> - Path to index file for input database.\n"+
                            "\t<search coverage percentage> - A number 0.0-1.0, represents how much of the database to search\n"+
                            "\t<query seqence> - The seqence to search for in the input db.");
        System.exit(1);
    }
}
