/**
 * This is a utility class that can represent an integer tuple
 * In the BLAST implementation, it is used to relate a word 
 * in a list with its index in the query
 * 
 * @author Sean McGroty
 */

public class HSP
{
	/**
	 * Basic constructor
	 * 
	 * @param q an integer
	 * @param s another integer
	 */
    public HSP(int q, int s)
    {
        qPos = q;
        sPos = s;
    }
    
    public int qPos;
    public int sPos;
}