public class BLASTN extends BLAST
{
    private int match;
    private int mismatch;
    
    public BLASTN()
    {
        //these are based on NCBI's defaults
        this.wordLength = 11;
        this.scoreCutoff = 25;
        this.match = 1
        this.mismatch = -3;
        this.gapOpenPenalty = 5;
        this.gapExtendPenalty = 2;
    }
    
    public BLASTN(int wordLength, int scoreCutoff, int match, int mismatch, int gapOpen, int gapExtend)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.match = match;
        this.mismatch = mismatch;
        this.gapOpenPenalty = gapOpen;
        this.gapExtendPenalty = gapExtend;
    }
    
    //Nucleotide words are just every possible word
    //analogous to protein version
    protected int[] findSeeds(String query)
    {
        int[] seeds = new int[query.length() - (wordLength-1)];
        for(int i = 0; i < seeds.length; i++)
            seeds[i] = i;
        
        return seeds;
    }
    
    protected int getScore(char a, char b)
    {
        if(a == b)
            return match;
        return mismatch;
    }
}