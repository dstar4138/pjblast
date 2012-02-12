/**
 * Parallel Version of BLAST.
 * Runs a search over a database list in parallel, searching each one on a separate process and reporting 
 * back the alignment ranges are. Also the parallel BLAST function makes use of a multicore environment to 
 * handle long sub word searches within the database entry itself.
 **/

import edu.rit.compbio.seq.Alignment;
import edu.rit.compbio.seq.ProteinDatabase; // Reads in FASTA style databases.
import edu.rit.compbio.seq.ProteinSequence;
import edu.rit.compbio.seq.Sequence;

import edu.rit.pj.Comm;
import edu.rit.pj.WorkerTeam;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerLongForLoop;
import edu.rit.mp.buf.ObjectItemBuf;

import java.io.File;

public class BlastRunnerPar {

    public static void main( String[] args ) throws Exception{
        Comm.init(args);
        final Comm world = Comm.world();
        final int rank = world.rank();
         
        if( args.length != 4 ){ usage(); }
        
        final ProteinDatabase pd = new ProteinDatabase( new File( args[1] ) , new File("tmpIndexFile.dat") );
        final float percentage = (new Float(args[2])).floatValue();
        
        final Sequence query;
        final BLAST aligner;
        if(args[0].toLowerCase().startsWith("p")){
            query = new ProteinSequence(">Search Query", args[3]);
            aligner = new BLASTP();
        }else{
            query = new NucleotideSequence(">Search Query", args[3]);
            aligner = new BLASTN();
        }
         
        // Loop over database and find alignments.
        new WorkerTeam().execute( new WorkerRegion(){
            public void run() throws Exception {
                execute( 0L,  (long) (pd.getDatabaseLength()*percentage), new WorkerLongForLoop(){
                    public void run( long first, long last ){
                        for( long i=first; i<=last; ++i){
                        	try{
	                            Alignment[] tmp = aligner.align( query, pd.getProteinSequence( i ) );
	                            if(rank != 0){
	                            	world.send(0, new ObjectItemBuf<Alignment[]>(tmp));
	                            }else{
	                            	//TODO: print alignments
	                            }
                        	}catch( Exception e ){}
                        }
                    }
               });
             }
        });
        
    }

    private static void usage(){
        System.out.println( "Usage: java BlastRunnerSeq [p/n] <input database> <search coverage percentage> <query seqence>\n"+
                            "\t[p/n] - p if proteins, n if nucleotides.\n"+
                            "\t<input database> - Path to a file with FASTA syle sequences.\n" +
                            "\t<search coverage percentage> - A number 1-100, represents how much of the database to search\n"+
                            "\t<query seqence> - The seqence to search for in the input db.");
        System.exit(1);
    }
}
