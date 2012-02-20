/**
 * Class BLAST is an abstract base class that implements the original BLAST ungapped
 * alignment algorithm. This class implements components common to both nucleotide and 
 * protein sequence alignments
 * 
 * @author: Sean McGroty
 */

import pjbio.Alignment;
import pjbio.Sequence;
import java.util.HashSet;
import java.util.ArrayList;

abstract public class BLAST
{
	protected int scoreCutoff;
	protected int wordLength;
	protected double eCutoff, K, LAM;
	protected long dbLength;

	/**
	 * Determines the score between a pair of residues
	 * 
	 * @param a the first residue
	 * @param b the second residue
	 * @return the score of a paired with b
	 */
	protected abstract int getScore(byte a, byte b);

	/**
	 * Finds the set of significant words from a query sequence
	 * 
	 * @param str a sequence to be decomposed into words
	 * @return an array of indices representing the position of each word
	 */
	protected abstract int[] findSeeds(byte[] str);

	/**
	 * Performs the BLAST ungapped alignment by finding all significant
	 * matches of the query in the subject
	 * 
	 * @param querySeq the query sequence
	 * @param subjectSeq the subject sequence
	 * @return An array of Alignments
	 */
	public Alignment[] align(Sequence querySeq, Sequence subjectSeq)
	{
		ArrayList<HSP> hits = new ArrayList<HSP>();
		HashSet<Alignment> alignments = new HashSet<Alignment>();
		int[] seeds;
		byte[] query;
		byte[] subject;
		int queryIndex, subjectIndex;
		double eScore;
		
		query = querySeq.sequence();
		subject = subjectSeq.sequence();
		
		//0: find the significant seeds
		seeds = findSeeds(query);

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
				pos = indexOf(subject,getRange(query,seeds[i],seeds[i] + wordLength-1),pos+1);
			}
			pos = 1;
		}   	

		//2: look for exact matches for seed words in the subject
		for(int i = 0; i < hits.size(); i++)
		{
			HSP temp = hits.get(i);
			int currScore, fScore = 0, rScore = 0, alignscore = 0, difference = 0, startRange = 0, endRange = 0;

			//extend forward
			queryIndex = temp.qPos + wordLength;
			subjectIndex = temp.sPos + wordLength;

			while(queryIndex < query.length && subjectIndex < subject.length)
			{
				fScore = alignscore;
				alignscore += getScore(query[queryIndex],subject[subjectIndex]);
				difference += (fScore - alignscore);

				//stop extending if the accumulated difference reaches the cutoff
				if(difference >= scoreCutoff) break;

				endRange++;
				queryIndex++;
				subjectIndex++;
			}

			//reset for next extension
			difference = 0;
			alignscore = 0;
			queryIndex = temp.qPos - 1;
			subjectIndex = temp.sPos - 1;

			//extend backwards
			while(queryIndex > 0 && subjectIndex > 0)
			{
				rScore = alignscore;
				alignscore += getScore(query[queryIndex],subject[subjectIndex]);
				difference += (rScore - alignscore);

				//stop extending if the accumulated difference reaches the cutoff
				if(difference >= scoreCutoff) break;

				startRange++;
				queryIndex--;
				subjectIndex--;
			}

			currScore = rawScore(getRange(query,temp.qPos - startRange,temp.qPos + wordLength + endRange),getRange(subject,temp.sPos - startRange,temp.sPos + wordLength + endRange));
			eScore = findEScore(currScore,query.length-1,dbLength);
			
			//3: keep only the extended alignments that pass cutoff 
			if((currScore > 0) && (eScore <= eCutoff))
			{	
				Alignment alignment = new Alignment();
				alignment.setMyQueryLength(query.length-1);
				alignment.setMySubjectLength(subject.length-1);
				alignment.setMyQueryStart(temp.qPos - startRange);
				alignment.setMyQueryFinish(temp.qPos + wordLength + endRange);
				alignment.setMySubjectStart(temp.sPos - startRange);
				alignment.setMySubjectFinish(temp.sPos + wordLength + endRange);
				alignment.setMyScore(currScore);
				alignments.add(alignment);
			}
		}

		return alignments.toArray(new Alignment[0]);
	}

	/**
	 * Computes the expect score for an alignment
	 * 
	 * @param score the alignment's raw score
	 * @param qLen the length of the aligned range in the query
	 * @param sLen the total length of the subject space (ex. the database size)
	 * @return the alignment's E score
	 */
	private double findEScore(int score, int qLen, long sLen)
	{
		double result;
		result = K * qLen * sLen * Math.exp(-score*LAM);	
		return result;
	}
	
	/**
	 * A basic search function that finds the first occurrence of pat in
	 * str, beginning at start
	 * 
	 * @param str a string
	 * @param pat the substring to be searched for
	 * @param start the start index in str
	 * @return the index of the first occurrence of pat in str, or -1 if not found
	 */
	private int indexOf(byte[] str, byte[] pat, int start)
	{
		int found = -1;

		for(int i = start; i < str.length - pat.length && found == -1; i++)
		{
			int j = 0;
			while(j < pat.length && str[i+j] == pat[j])
				j++;
			if(j == pat.length)
				found = i;
		}
		return found;
	}

	/**
	 * Returns a subrange of a byte array, from start to end inclusive
	 * 
	 * @param array input array
	 * @param from start index
	 * @param to end index
	 * @return a subarray of array as indicated by from and to
	 */
	private byte[] getRange( byte[] array, int from, int to){
		byte[] newarray = new byte[ (to - from)];
		for(int i=from,j=0; i<to; i++,j++) newarray[j]=array[i];
		return newarray;
	} 

	/**
	 * Calculates the raw score (sum of all pairwise scores)
	 * of two aligned strings, as scored by the scoring function
	 * 
	 * @param a first string
	 * @param b second string
	 * @return the raw score
	 */
	int rawScore(byte[] a, byte[] b)
	{
		int score = 0;
		for(int i = 0; i < a.length && i < b.length; i++)
		{
			score += getScore(a[i],b[i]);
		}
		return score;
	}
}

