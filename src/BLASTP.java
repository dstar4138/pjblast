/**
 * Class BLASTP implements protein sequence alignments. It provides the protein-specific
 * word generation and scoring functions used by class BLAST.
 * 
 * @author Sean McGroty
 */

import java.util.ArrayList;
import pjbio.Blosum62;

public class BLASTP extends BLAST
{    
    private int[][] scoringMatrix;
    private int wordCutoff;
    
    /**
     * The default constructor sets parameters for use with the BLOSUM62 matrix
     * @param dbLength the total number of amino acids in the database
     */
    public BLASTP(long dbLength)
    {
        //these are based on NCBI's defaults
        this.wordLength = 3;
        this.wordCutoff = 13;
        this.scoreCutoff = 23;
        this.scoringMatrix = Blosum62.matrix;
        this.eCutoff = 10;
        this.K = 0.134;	//values obtained from BLAST, page 66
        this.LAM = 0.318;
        this.dbLength = dbLength;
    }
    
    /**
     * Constructor for using a user-provided scoring matrix. All parameters must be customized for that matrix
     * since the use of another scoring matrix strongly implies differing statistical measures for scoring and analysis
     * 
     * @param wordLength
     * @param wordCutoff
     * @param scoreCutoff
     * @param eCut
     * @param userScoringMatrix
     */
    public BLASTP(int wordLength, int wordCutoff, int scoreCutoff, double eCut, int[][] userScoringMatrix)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.scoringMatrix = userScoringMatrix;
        this.wordCutoff = wordCutoff;
        this.eCutoff = eCut;
    }
    
    /**
     * getScore applies the scoring matrix to the two amino acids represented
     * by a and b
     * 
     * @param a first amino acid
     * @param b second amino acid
     * @return the score from the matrix
     */
    protected int getScore(byte a, byte b)
    {
        return scoringMatrix[a][b];
    }
    
    /**
     * Finds high scoring words by scanning each word along the query
     * and returning only ones that yield a pairwise score above a cutoff
     * 
     * @param query the query 
     * @return an integer array of indexes into the query representing high scoring words
     */
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
      
		int[] retval = new int[foundSeeds.size()];
		
		for(int i = 0; i < retval.length; i++)
			retval[i] = foundSeeds.get(i).intValue();
		return retval;
    }   
}