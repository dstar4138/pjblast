import java.util.ArrayList;

import modpj.Blosum62;


public class BLASTP extends BLAST
{    
    private int[][] scoringMatrix;
    
    private int wordCutoff;
    
    public BLASTP()
    {
        //these are based on NCBI's defaults
        this.wordLength = 3;
        this.wordCutoff = 13;
        this.scoreCutoff = 23;
        this.gapOpenPenalty = 11;
        this.gapExtensionPenalty = 1;
        this.scoringMatrix = Blosum62.matrix;
        this.eCutoff = 10;
        this.K = 0.035;
        this.LAM = 0.252;
    }
    
    //TODO: add invalid scoring matrix error handling
    public BLASTP(int wordLength, int wordCutoff, int scoreCutoff,int gapOpen, int gapExtend, double eCut, int[][] userScoringMatrix)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.gapOpenPenalty = gapOpen;
        this.gapExtensionPenalty = gapExtend;
        this.scoringMatrix = userScoringMatrix;
        this.wordCutoff = wordCutoff;
        this.eCutoff = eCut;
    }
    
    protected int getScore(byte a, byte b)
    {
        return scoringMatrix[a][b];
    }
    
    //self-scans the query for high-scoring words
    //NOTE: this is currently not a very efficient way of doing this and is O(n^2*wordLength)
    protected int[] findSeeds(byte[] query)
    {
        ArrayList<Integer> foundSeeds = new ArrayList<Integer>();
        int currScore;
        
        //this loop walks through the query, breaking it into words
        for(int i = 1; i <= query.length - wordLength; i++)
        {
            currScore = 0;
            //this loops steps through the full length of the query
            for(int j = 1; j <= query.length - wordLength; j++)
            {
                currScore = 0;
                //score the current word across the entire query
                for(int k = 0; k < wordLength; k++)
                {
                    currScore += getScore(query[i+k],query[j+k]);
                }
                
                //if the word obtains a sufficient score against any portion of the query, keep it and stop
                //i: the high-scoring word as indicated by index
                if(currScore >= wordCutoff)
                {
                    foundSeeds.add(new Integer(i));
                    break;
                }
            }
        }
        //create Alignment objects
        
        
		int[] retval = new int[foundSeeds.size()];
		
		for(int i = 0; i < retval.length; i++)
			retval[i] = foundSeeds.get(i).intValue();
		return retval;
    }
    
  /*  protected Alignment[] doGapped(AlignRange[] range, Sequence query, Sequence subject)
    {
    	SequencePair[] ungapped = new SequencePair[range.length];
    	
        for(int i = 0; i < ungapped.length; i++)
        {
        	//System.out.println("qrange "+ (range[i].qStart-1) + "," + range[i].qEnd + "srange "+ (range[i].sStart-1) + "," + range[i].sEnd);
        	ungapped[i] = new SequencePair(new ProteinSequence(query.description(),query.elementsToString().substring(range[i].qStart-1, range[i].qEnd)),
        								   new ProteinSequence(subject.description(),subject.elementsToString().substring(range[i].sStart-1, range[i].sEnd)));
        }
    	
        
        Alignment[] results = new Alignment[ungapped.length];
    	ProteinLocalAlignmentSeq aligner = new ProteinLocalAlignmentSeq();
    	//aligner stuff
    	
    	for(int i = 0; i < ungapped.length; i++)
        {
        	//System.out.println(ungapped[i].query.elementsToString());
        	//System.out.println(ungapped[i].subject.elementsToString());
    		aligner.setQuerySequence((ProteinSequence) ungapped[i].query, (long)i);
    		aligner.setSubjectSequence((ProteinSequence) ungapped[i].subject, (long)i);
            results[i] = aligner.align();
        }
    	
    	return results;
    }*/
}