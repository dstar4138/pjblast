/**
 * Sequential Version of BLAST.
 * Only run on smaller databases of sequences, as the sequential version could theoretically 
 * take exponential time. 
 **/

import edu.rit.compbio.seq.ProteinDatabase; // Reads in FASTA style databases.

import java.io.File;

public class BlastRunnerSeq {

    public static void main( String[] args ){
        
        if( args.length != 3 ){ usage(); }

        String inputDatabase = args[0];
        float percentage = Float.floatValue(args[1]);
        String query = args[2];
        
        ProteinDatabase pd = ProteinDatabase( new File( inputDatabase ) , new File("tmpIndexFile.dat") );
        BLASTP aligner = new BLASTP();

        for(long i=0; i < pd.getDatabaseLength()*percentage; i++){
            AlignRange[] tmp = aligner.align( query, pd.getProteinSequence( i ).toString() );
           // print the alignment nicely.  
        }
    }

    private static void usage(){
        System.out.println( "Usage: java BlastRunnerSeq <input database> <search coverage percentage> <query seqence>\n"+
                            "\t<input database> - Path to a file with FASTA syle sequences.\n" +
                            "\t<search coverage percentage> - A number 1-100, represents how much of the database to search\n"+
                            "\t<query seqence> - The seqence to search for in the input db.");
        System.exit(1);
    }

}
