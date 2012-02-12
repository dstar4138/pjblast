import edu.rit.compbio.seq.Sequence;

public class SequencePair {
	public SequencePair(Sequence q, Sequence s)
	{
		query = q;
		subject = s;
	}
	
	public Sequence query;
	public Sequence subject;
}
