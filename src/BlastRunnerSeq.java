/**
 * Sequential Version of BLAST.
 * Only run on smaller databases of sequences, as the sequential version could theoretically 
 * take exponential time. 
 **/

import edu.rit.compbio.seq.ProteinDatabase; // Reads in FASTA style databases.
import edu.rit.compbio.seq.ProteinSequence;
import edu.rit.compbio.seq.Sequence;
import edu.rit.pj.Comm;

import java.io.File;

public class BlastRunnerSeq {

    public static void main( String[] args ){
        Comm.init(args);
         
        if( args.length != 4 ){ usage(); }
        
        ProteinDatabase pd = new ProteinDatabase( new File( args[1] ) , new File("tmpIndexFile.dat") );
        float percentage = (new Float(args[2])).floatValue();
        
        Sequence query;
        BLAST aligner;
        if(args[0].toLowerCase().startsWith("p")){
            query = new ProteinSequence("Search Query", args[3]);
            aligner = new BLASTP();
        }else{
            query = new NucleotideSequence("Search Query", args[3]);
            aligner = new BLASTN();
        }
        
        for(long i=0; i < pd.getDatabaseLength()*percentage; i++){
            AlignRange[] tmp = aligner.align( query, pd.getProteinSequence( i ) );
            if(tmp.length>0){
                //TODO: print the alignment nicely, for this subject sequence vs the query.
            }
        }
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
