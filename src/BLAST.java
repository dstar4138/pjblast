
import java.util.ArrayList;
import java.util.Arrays;
import modpj.Alignment;
import modpj.Sequence;

abstract public class BLAST
{
    protected int scoreCutoff;
    protected int wordLength;
    protected int gapOpenPenalty;
    protected int gapExtensionPenalty;
    protected double eCutoff, K, LAM;
    protected String queryDesc, subjectDesc;
    
    //returns the score for any two letter pairs
    protected abstract int getScore(byte a, byte b);
    
    //creates the list of words to be used in the seeding
    protected abstract int[] findSeeds(byte[] str);
       
    public Alignment[] align(Sequence querySeq, Sequence subjectSeq)
    {
    	queryDesc = querySeq.description();
    	subjectDesc = subjectSeq.description();
        byte[] query = querySeq.sequence();
        byte[] subject = subjectSeq.sequence();
        //seeds is an array of indexes into the query representing the words
        int[] seeds = findSeeds(query);
        ArrayList<HSP> hits = new ArrayList<HSP>();
		ArrayList<Alignment> alignments = new ArrayList<Alignment>();
        int alignscore = 0;
        //double eScore;
        int startRange, endRange, queryIndex, subjectIndex;
        
        //1: find all exact matches between a word and some position in the subject
        //start from 1 position, 0 position is unused
        int pos = 1;
        for(int i = 0; i < seeds.length; i++)
        {
        	pos = indexOf(subject,getRange(query,seeds[i],seeds[i] + wordLength),pos);
        	
            while(pos != -1)
            {
                //HSP(word's position in query, word's position in subject)
                hits.add(new HSP(seeds[i],pos));
                pos = indexOf(subject,getRange(query,seeds[i],seeds[i] + wordLength),pos+1);
            }
			pos = 1;
        }
        
        //2: extend the match forwards and backwards until the score decreases too much
        for(int i = 0; i < hits.size(); i++)
        {
			HSP temp = hits.get(i);
            int currScore = 0;
            int difference = 0;
            endRange = 0;
            
            //extend forward
            queryIndex = temp.qPos + wordLength;
            subjectIndex = temp.sPos + wordLength;
            
            while(queryIndex < query.length && subjectIndex < subject.length)
            {
                currScore = alignscore;
                alignscore += getScore(query[queryIndex],subject[subjectIndex]);
                difference += (currScore - alignscore);
                
                //stop extending if the accumulated difference reaches the cutoff
                if(difference >= scoreCutoff) break;
                
                endRange++;
                queryIndex++;
                subjectIndex++;
            }
            
            //reset for next extension
            currScore = 0;
            difference = 0;
            queryIndex = temp.qPos - 1;
            subjectIndex = temp.sPos - 1;
            startRange = 0;
            
            //extend backwards
            while(queryIndex > 0 && subjectIndex > 0)
            {
                currScore = alignscore;
                alignscore += getScore(query[queryIndex],subject[subjectIndex]);
                difference += (currScore - alignscore);
                
                //stop extending if the accumulated difference reaches the cutoff
                if(difference >= scoreCutoff) break;
                
                startRange++;
                queryIndex--;
                subjectIndex--;
            }
            
            //add check for identical ranges
			
			//3: keep only the extended alignments that pass cutoff 
            System.out.println(findEScore(currScore,querySeq.length(),subjectSeq.length()));
            if(findEScore(currScore,querySeq.length(),subjectSeq.length()) <= eCutoff)
            {
                Alignment alignment = new Alignment();
                alignment.setMyQueryId(0L);
            	alignment.setMySubjectId(0L);
            	alignment.setMyQueryLength(querySeq.length());
            	alignment.setMySubjectLength(subjectSeq.length());
            	alignment.setMyQueryStart(temp.qPos - startRange);
            	alignment.setMyQueryFinish(temp.qPos + (wordLength-1) + endRange);
            	alignment.setMySubjectStart(temp.sPos - startRange);
            	alignment.setMySubjectFinish(temp.sPos + (wordLength-1) + endRange);
            	alignment.setMyScore(currScore);
            	alignment.setMyTraceback(new byte[alignment.getQueryFinish() - alignment.getQueryStart()]);
            	Arrays.fill(alignment.getMyTraceback(), (byte) Alignment.QUERY_ALIGNED_WITH_SUBJECT);    
                alignments.add(alignment);
            }
        }
    
        return alignments.toArray(new Alignment[0]);
    }
    
    private double findEScore(int score, int qLen, int sLen)
    {
        double y = K * qLen * sLen * Math.pow(Math.E,-(LAM*score));
        
        return 1 - Math.pow(Math.E,-y);
    }
    
    //naive search for now, could be improved with a DFA or some type of hash based search
    private int indexOf(byte[] str, byte[] pat, int start)
    {
    	int found = -1;
    	for(int i = start; i < str.length - pat.length; i++)
    	{
    		for(int j = 0; j < pat.length; j++)
    		{
    			if(str[i] != pat[j])
    				break;
    			found = i;
    		}
    	}
    	return found;
    }
    
    private byte[] getRange( byte[] array, int from, int to){
    	byte[] newarray = new byte[ (to - from)];
    	for(int i=from,j=0; i<to; i++,j++) newarray[j]=array[i];
    	return newarray;
    }                        
}
		
		