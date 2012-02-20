/**
 * Sequential Version of BLAST - 
 *
 * @author Alexander Dean
 **/

import edu.rit.pj.Comm;
import java.io.File;

import pjbio.*;

public class BlastRunnerSeq {

    public static void main( String[] args ) throws Exception{
        Comm.init(args);
         
        if( args.length != 5 ){ usage(); }

        AlignmentPrinter out;
        ProteinDatabase pd = new ProteinDatabase( new File( args[1] ) , new File(args[2]) );
        float percentage = (new Float(args[3])).floatValue();
        
        Sequence query;
        BLAST aligner;
        if(args[0].toLowerCase().startsWith("p")){
            query = new ProteinSequence(">Search Query", args[4]);
            aligner = new BLASTP(pd.getProteinCount());
        }else{
            query = new NucleotideSequence(">Search Query", args[4]);
            aligner = new BLASTN(pd.getProteinCount());
        }
        
        long t1 = System.currentTimeMillis();
        for(long i=0; i < pd.getProteinCount()*percentage; i++){
        	Alignment[] tmp = aligner.align( query, pd.getProteinSequence( i ) );
        	// Print the alignments
            for(int j=0; j<tmp.length;j++){
            	 out = new AlignmentPrinter(System.out,new DefaultAlignmentStats(tmp[j].getSubjectLength()));
            	 out.printDetails(tmp[j], query, pd.getProteinSequence(i));
            }
        }
        long t2 = System.currentTimeMillis();
        
        System.out.printf("Running Time: %d msec%n", t2-t1);
    }

    private static void usage(){
        System.out.println( "Usage: java BlastRunnerSeq [p/n] <input database> <index file> <search coverage percentage> <query seqence>\n"+
                            "\t[p/n] - p if proteins, n if nucleotides.\n"+
                            "\t<database indexs> - Path to index file for input database.\n"+
                            "\t<input database> - Path to a file with FASTA style sequences.\n" +
                            "\t<search coverage percentage> - A number 0.0-1.0, represents how much of the database to search\n"+
                            "\t<query seqence> - The sequence to search for in the input db.");
        System.exit(1);
    }

}
