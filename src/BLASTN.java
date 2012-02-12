import edu.rit.compbio.seq.Alignment;
import edu.rit.compbio.seq.Sequence;

public class BLASTN extends BLAST
{
    private int match;
    private int mismatch;
    
    public BLASTN()
    {
        //these are based on NCBI's defaults
        this.wordLength = 11;
        this.scoreCutoff = 25;
        this.match = 1;
        this.mismatch = -3;
        this.gapOpenPenalty = 5;
        this.gapExtensionPenalty = 2;
        this.eCutoff = 10;
        this.K = 0.13;
        this.LAM = 0.318;
    }
    
    public BLASTN(int wordLength, int scoreCutoff,int gapOpen, int gapExtend, double eCut)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.gapOpenPenalty = gapOpen;
        this.gapExtensionPenalty = gapExtend;
        this.eCutoff = eCut;
    }
    
    //Nucleotide words are just every possible word
    //output analogous to protein version
    protected int[] findSeeds(byte[] query)
    {
        int[] seeds = new int[query.length - wordLength];
        for(int i = 1; i < seeds.length; i++)
            seeds[i] = i;
        
        return seeds;
    }
    
    protected int getScore(byte a, byte b)
    {
        if(a == b)
            return match;
        return mismatch;
    }

	protected Alignment[] doGapped(AlignRange[] range, Sequence query,
			Sequence subject) {
		// TODO Auto-generated method stub
		return null;
	}
}