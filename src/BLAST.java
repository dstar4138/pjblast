import java.util.ArrayList;

abstract public class BLAST
{
    protected int scoreCutoff;
    protected int wordLength;
    protected int gapOpenPenalty;
    protected int gapExtensionPenalty;
    protected double eCutoff, K, LAM;
    
    //returns the score for any two letter pairs
    protected abstract int getScore(char a, char b);
    
    //creates the list of words to be used in the initial ungapped alignment
    protected abstract int[] findSeeds(String query);
    
    public AlignRange[] align(String query, String subject)
    {
        //seeds is an array of indexes into the query representing the words
        int[] seeds = findSeeds(query);
        ArrayList<HSP> hits = new ArrayList<HSP>();
		ArrayList<AlignRange> alignments = new ArrayList<AlignRange>();
        int alignscore = 0;
        double eScore;
        
        int startRange, endRange, queryIndex, subjectIndex;
        
        //1: find all exact matches between a word and some position in the subject
        int pos = 0;
        for(int i = 0; i < seeds.length; i++)
        {
            pos = subject.indexOf(query.substring(seeds[i],seeds[i] + wordLength),pos);
            while(pos != -1)
            {
                //HSP(word's position in query, word's position in subject)
                hits.add(new HSP(seeds[i],pos));
                pos = subject.indexOf(query.substring(seeds[i],seeds[i] + wordLength),pos+1);
            }
			pos = 0;
        }
        
        //2: extend the match forwards and backwards until the score decreases too much
        for(int i = 0; i < hits.size(); i++)
        {
			HSP temp = hits.get(i);
            int currScore = 0;
            endRange = 0;
            //extend forward
            queryIndex = temp.qPos + wordLength;
            subjectIndex = temp.sPos + wordLength;
            
            while(queryIndex < query.length() && subjectIndex < subject.length())
            {
                currScore = alignscore;
                alignscore += getScore(query.charAt(queryIndex),subject.charAt(subjectIndex));
                
                //stop extending if we get a new score that is lower than previous-cutoff
                if(alignscore < currScore - scoreCutoff) break;
                
                endRange++;
                queryIndex++;
                subjectIndex++;
            }
            
            //reset for next extension
            currScore = 0;
            queryIndex = temp.qPos - 1;
            subjectIndex = temp.sPos - 1;
            startRange = 0;
            
            //extend backwards
            while(queryIndex > -1 && subjectIndex > -1)
            {
                currScore = alignscore;
                alignscore += getScore(query.charAt(queryIndex),subject.charAt(subjectIndex));
                
                //stop extending if we get a new score that is lower than previous-cutoff
                if(alignscore < currScore - scoreCutoff) break;
                
                startRange++;
                queryIndex--;
                subjectIndex--;
            }
			
			//3: keep only the extended alignments that pass cutoff
                //query start, query end, subject start, subject end
                //compute the E score and keep only if above e score cutoff
            System.out.println("qrange " + (temp.qPos - startRange) + "," + (temp.qPos + (wordLength-1) + endRange));
            System.out.println("srange " + (temp.sPos - startRange) + "," + (temp.sPos + (wordLength-1) + endRange));
           
            eScore = findEScore(currScore, query.length(), subject.length());
            if(eScore >= eCutoff)
                alignments.add(new AlignRange(eScore, temp.qPos - startRange,temp.qPos + (wordLength-1) + endRange,temp.sPos - startRange,temp.sPos + (wordLength-1) + endRange));
        }
        
        AlignRange[] results = new AlignRange[alignments.size()];
        alignments.toArray(results);
        return results;
    }
    
    private double findEScore(int score, int qLen, int sLen)
    {
        double y = K * qLen * sLen * Math.pow(Math.E,-(LAM*score));
        
        return 1 - Math.pow(Math.E,-y);
    }
                             
}
		
		
