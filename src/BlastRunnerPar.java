/**
 * Parallel Version of BLAST.
 * Runs a search over a database list in parallel, searching each one on a separate process and reporting 
 * back the alignment ranges are. Also the parallel BLAST function makes use of a multicore environment to 
 * handle long sub word searches within the database entry itself.
 **/

import edu.rit.pj.Comm;
import edu.rit.pj.WorkerTeam;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerLongForLoop;

import java.io.File;

import modpj.*;


public class BlastRunnerPar {

	static ProteinDatabase pd;
	
    public static void main( String[] args ) throws Exception{
        Comm.init(args);
        //final Comm world = Comm.world();
        //final int rank = world.rank();
        //final int size = world.size();
         
        if( args.length != 5 ){ usage(); }
        // TODO: ProteinDatabase wont necessarily work if the user passes in a nucleotide file in.
        pd = new ProteinDatabase( new File( args[1] ) , new File(args[2]) );
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
        new WorkerTeam().execute( new WorkerRegion(){
            public void run() throws Exception {
                execute( 0L,  (long) (pd.getProteinCount()*percentage), new WorkerLongForLoop(){
                    public void run( long first, long last ){
                    	AlignmentPrinter out;
                        for( long i=first; i<=last; ++i){
                        	try{
	                            Alignment[] tmp = aligner.align( query, pd.getProteinSequence( i ) );
	                            
	                            out = new AlignmentPrinter(System.out,new DefaultAlignmentStats(pd.getProteinSequence(i).length()));
	                    		for(Alignment align : tmp)
	                       	 		out.printDetails(align, query, pd.getProteinSequence( i ));
	                           
                        	}catch( Exception e ){
                        		//LATER: Logging?
                        	}
                        }
                    }
               });
             }
        });
        
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
