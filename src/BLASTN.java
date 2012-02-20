/**
 * Class BLASTN implements nucleotide alignments. It provides the scoring and word generation functions
 * used in nucleotide sequence alignments
 * @author Sean McGroty
 */

public class BLASTN extends BLAST
{
    private int match;
    private int mismatch;
    
    /**
     * Basic constructor for nucleotide BLAST
     * 
     * @param dbLength total number of nucleotides in the search space
     * 
     */
    public BLASTN(long dbLength)
    {
        //these are based on NCBI's defaults
        this.wordLength = 11;
        this.scoreCutoff = 25;
        this.match = 1;
        this.mismatch = -2;
        this.eCutoff = 10;
        this.K = 0.182;
        this.LAM = 0.192;	//WU-BLAST default scores for ungapped nucleotide alignment
        this.dbLength = dbLength;
    }
    
    /**
     * Constructor allowing user-specified parameters
     * 
     * @param wordLength
     * @param scoreCutoff
     * @param match
     * @param mismatch
     * @param eCut
     * @param dbLength
     * @param K
     * @param lambda
     */
    public BLASTN(int wordLength, int scoreCutoff, int match, int mismatch, double eCut, int dbLength, double K, double lambda)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.match = match;
        this.mismatch = mismatch;
        this.eCutoff = eCut;
        this.K = K;
        this.LAM = lambda;
        this.dbLength = dbLength;
    }
    
    /**
     * Returns nucleotide seeds, which is just every possible word
     *
     * @param query Nucleotide query sequence
     * @return array consisting of indexes for all possible words
     */
    protected int[] findSeeds(byte[] query)
    {
        int[] seeds = new int[query.length - wordLength];
        for(int i = 1; i < seeds.length; i++)
            seeds[i] = i;
        
        return seeds;
    }
    
    /**
     * Returns match/mismatch scores
     * 
     * @param a first nucleotide
     * @param b second nucleotide
     * @return score corresponding to match or mismatch
     */
    protected int getScore(byte a, byte b)
    {
        if(a == b)
            return match;
        return mismatch;
    }
}