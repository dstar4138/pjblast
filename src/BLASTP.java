public class BLASTP extends BLAST
{
    private final int[][] BLOSUM62 = 
       {{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
        {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
        {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
        {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
        { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
        {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
        {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
        { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
        {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
        {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
        {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
        {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
        {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
        {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
        { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
        { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
        {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
        { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};
    
    private int[][] scoringMatrix;
    
    public BLASTP()
    {
        //these are based on NCBI's defaults
        this.wordLength = 3;
        this.scoreCutoff = 25;
        this.gapOpenPenalty = 11;
        this.gapExtendPenalty = 1;
        this.scoringMatrix = BLOSUM62;
    }
    
    //TODO: add invalid scoring matrix error handling
    public BLASTP(int wordLength, int scoreCutoff,int gapOpen, int gapExtend, int[][] userScoringMatrix)
    {
        this.wordLength = wordLength;
        this.scoreCutoff = scoreCutoff;
        this.gapOpenPenalty = gapOpen;
        this.gapExtendPenalty = gapExtend;
        this.scoringMatrix = userScoringMatrix;
    }
    
    private int getScore(char a, char b)
    {
        return scoringMatrix[aaToN(res1)][aaToN(res2)];
    }
    
    //self-scans the query for high-scoring words
    //NOTE: this is currently not a very efficient way of doing this and is O(n^2*wordLength)
    private int[] findSeeds(String query)
    {
        ArrayList foundSeeds = new ArrayList();
        int currScore = 0;
        
        //this loop walks through the query, breaking it into words
        for(int i = 0; i < query.length() - wordLength; i++)
        {
            //this loops steps through the full length of the query
            for(int j = 0; j < query.length() - wordLength, j++)
            {
                //score the current word across the entire query
                for(int k = 0; k < wordLength; k++)
                {
                    currScore = getScore(query.charAt(i+k),subject.charAt(j+k));
                }
                //i: the high-scoring word as indicated by index
                if(currScore >= scoreCutoff)
                    foundSeeds.add(new Integer(i));
            }
        }
        return foundSeeds.toArray();
    }    
    
    //TODO: invalid character error handling
    private int aaToN(char c) {
        switch(Character.toLowerCase(c))
        {
            case 'a':
                return 0;
            case 'r':
                return 1;
            case 'n':
                return 2; 
            case 'd':
                return 3;  
            case 'c':
                return 4;
            case 'q':
                return 5; 
            case 'e':
                return 6;   
            case 'g':
                return 7;
            case 'h':
                return 8;
            case 'i':
                return 9;
            case 'l':
                return 10;  
            case 'k':
                return 11;     
            case 'm':
                return 12;   
            case 'f':
                return 13;    
            case 'p':
                return 14;    
            case 's':
                return 15;
            case 't':
                return 16;    
            case 'w':
                return 17;
            case 'y':
                return 18;
            case 'v':
                return 19;
        }
        return -1;
    }
}