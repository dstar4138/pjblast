public class BLAST
{
    protected int scoreCutoff;
    protected int wordLength;
    protected int gapOpenPenalty;
    protected int gapExtensionPenalty;
    
    //returns the score for any two letter pairs
    protected abstract int getScore(char a, char b);
    
    //creates the list of words to be used in the initial ungapped alignment
    protected abstract int[] findSeeds(String query);
    
    protected void align(String subject, String query)
    {
        //seeds is an array of indexes into the query representing the words
        int[] seeds = findSeeds(subject, query);
        ArrayList hits = new ArrayList();
        
        int startRange, endRange;
        
        //1: find all exact matches between a word and some position in the query
        int pos;
        for(int i = 0; i < seeds.length; i++)
        {
            pos = subject.indexOf(query.substring(seeds[i],seeds[i] + (wordLength-1)),pos);
            while(pos != -1)
            {
                //HSP(matching word index, position in query)
                hits.add(new HSP(i,pos));
                pos = subject.indexOf(query.substring(seeds[i],seeds[i] + (wordLength-1)),pos);
            }
        }
        
        //2: extend the match forwards and backwards until the score decreases
        for(int i = 0; i < hits.length(); i++)
        {
            //extend forward
            for(endRange = 0; endRange < subject.length() - wordLength && j < query.length() - wordLength; endRange++)
            {
                //stop extending if we get a negative score
                if(getScore(query.charAt(endRange+(hits.get(i).qPos + wordLength),subject.charAt(endRange+(hits.get(i).sPos + wordLength)))) < 0)
                    break;
            }
            
            //extend backwards
            for(startRange = 0; startRange < 
        }

        //3: keep only the extended alignments that pass E score
        //4: still working
    }
                             
}
		
		