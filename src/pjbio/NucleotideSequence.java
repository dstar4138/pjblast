package pjbio;

/**
 * Since we are also implementing BLASTN, we are providing a nucleotide sequence 
 * since the ProteinSequence was already implemented by the ParallelJava library.
 **/

public class NucleotideSequence extends Sequence{

	private String savedSequence;
	
	public NucleotideSequence( String description, String sequence ){
    	super();
    	myDescription = description;
    	savedSequence = sequence;
    }

	@Override
	public char charAt(int arg0) {
		return savedSequence.charAt(arg0);
	}

	@Override
	public String toString(){
		return "NucleotideSequence( " + myDescription + " )";
	}

}
