import modpj.Alignment;
import modpj.AlignmentPrinter;
import modpj.DefaultAlignmentStats;
import modpj.ProteinSequence;


public class ModuleTest extends BLASTP
{
    public static void main(String[] args)
    {
        BLASTP aligner = new BLASTP();
        ProteinSequence q, s;
        
        switch(Integer.parseInt(args[0]))
        {
            case 1:
            	q = new ProteinSequence(">Query", args[1]);
                int[] seeds = aligner.findSeeds(q.sequence());
                
                System.out.println("Seed indexes for " + q.description());
                
                for(int i = 0; i < seeds.length; i++)
                {
                    System.out.print(seeds[i] + " ");
                }
                System.out.println();
                break;
                
            case 2:
            	q = new ProteinSequence(">Query",args[1]);
            	s = new ProteinSequence(">Subject",args[2]);
            	Alignment[] results = aligner.align(q ,s);
                AlignmentPrinter out = new AlignmentPrinter(System.out,new DefaultAlignmentStats(s.length()));
                for(int i = 0; i < results.length; i++)
                {
                    out.printDetails(results[i], q, s);
                }
                break;
            case 3:
              /*  int score = 0;
                for(int i = 0; i < args[1].length(); i++)
                {
                    System.out.print(aligner.getScore(args[1].charAt(i),args[2].charAt(i)) + " ");
                    score += aligner.getScore(args[1].charAt(i),args[2].charAt(i));
                }
                System.out.println("Total " + score);*/
            	break;
            case 4:
            	/*q = new ProteinSequence(">Query",args[1]);
            	s = new ProteinSequence(">Subject",args[2]);
            	ProteinLocalAlignmentSeq test = new ProteinLocalAlignmentSeq();
                test.setQuerySequence(q,1L);
                test.setSubjectSequence(s,1L);
                Alignment output = test.align();
                AlignmentPrinter out1 = new AlignmentPrinter(System.out,new DefaultAlignmentStats(s.length()));
                out1.printDetails(output, q, s);*/
            	break;
        }
        
    }
    
    public static void usage()
    {
        System.out.println("Usage: java ModuleTest testType query subject");
    }
}