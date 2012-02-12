//Represents two corresponding ranges and their score

public class AlignRange
{
    public AlignRange(double score, int q1, int q2, int s1, int s2)
    {
        qStart = q1;
		qEnd = q2;
		sStart = s1;
		sEnd = s2;
        eScore = score;
    }
    
    public int qStart, qEnd, sStart, sEnd;
    public double eScore;
}