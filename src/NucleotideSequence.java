/**
 * Since we are also implementing BLASTN, we are providing a nucleotide sequence 
 * since the ProteinSequence was already implemented by the ParallelJava library.
 **/

import edu.rit.compbio.seq.Sequence;
public class NucleotideSequence extends Sequence{

    private String myDescription;
	private String mySequence;

	public NucleotideSequence( String description, String sequence ){
    	super();
    	myDescription = description;
    	mySequence = sequence;
    }

	@Override
	public char charAt(int arg0) {
		return mySequence.charAt(arg0);
	}

	@Override
	public String toString(){
		return "NucleotideSequence( " + myDescription + " )";
	}

}
