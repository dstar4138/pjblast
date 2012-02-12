import edu.rit.compbio.seq.Sequence;
import edu.rit.compbio.seq.Alignment;
import java.util.ArrayList;
import java.util.Arrays;

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
    
    //creates the list of words to be used in the initial ungapped alignment
    protected abstract int[] findSeeds(byte[] str);
    
    protected abstract Alignment[] doGapped(AlignRange[] range, Sequence query, Sequence subject);
    
    public Alignment[] align(Sequence querySeq, Sequence subjectSeq)
    {
    	queryDesc = querySeq.description();
    	subjectDesc = subjectSeq.description();
        byte[] query = querySeq.sequence();
        byte[] subject = subjectSeq.sequence();
        //seeds is an array of indexes into the query representing the words
        int[] seeds = findSeeds(query);
        ArrayList<HSP> hits = new ArrayList<HSP>();
		ArrayList<AlignRange> alignments = new ArrayList<AlignRange>();
        int alignscore = 0;
        double eScore;
        int startRange, endRange, queryIndex, subjectIndex;
        
        //1: find all exact matches between a word and some position in the subject
        //start from 1 position, 0 position is unused
        int pos = 1;
        for(int i = 0; i < seeds.length; i++)
        {
        	pos = indexOf(subject,Arrays.copyOfRange(query,seeds[i],seeds[i] + wordLength),pos);
        	
            while(pos != -1)
            {
                //HSP(word's position in query, word's position in subject)
                hits.add(new HSP(seeds[i],pos));
                pos = indexOf(subject,Arrays.copyOfRange(query,seeds[i],seeds[i] + wordLength),pos+1);
            }
			pos = 1;
        }
        
        //Exact match positions
        for(int i = 0; i < hits.size(); i++)
        {
        	System.out.println(hits.get(i).qPos + "," + hits.get(i).sPos);
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
			
			//3: keep only the extended alignments that pass cutoff            
            //debug statements
            System.out.println("qrange " + (temp.qPos - startRange) + "," + (temp.qPos + (wordLength-1) + endRange));
            System.out.println("srange " + (temp.sPos - startRange) + "," + (temp.sPos + (wordLength-1) + endRange));
           
            //-1 to length to account for unused 0 position
            eScore = findEScore(currScore, query.length - 1, subject.length - 1);
            //if(eScore >= eCutoff)
                //E score, query start, query end, subject start, subject end
                alignments.add(new AlignRange(eScore, temp.qPos - startRange,temp.qPos + (wordLength-1) + endRange,temp.sPos - startRange,temp.sPos + (wordLength-1) + endRange));
        }
        
        //AlignRange[] temp = alignments.toArray(temp);
        Alignment[] results = doGapped((AlignRange[]) alignments.toArray(new AlignRange[0]),querySeq,subjectSeq);
    
        return results;
    }
    
    private double findEScore(int score, int qLen, int sLen)
    {
        double y = K * qLen * sLen * Math.pow(Math.E,-(LAM*score));
        
        return 1 - Math.pow(Math.E,-y);
    }
    
    //naive search for now
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
                             
}
		
		