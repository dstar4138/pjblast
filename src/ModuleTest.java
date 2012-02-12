public class ModuleTest extends BLASTP
{
    public static void main(String[] args)
    {
        BLASTP aligner = new BLASTP();
        
        switch(Integer.parseInt(args[0]))
        {
            case 1:
                int[] seeds = aligner.findSeeds(args[1]);
                System.out.println("Seed indexes for seq1");
                
                for(int i = 0; i < seeds.length; i++)
                {
                    System.out.print(seeds[i] + " ");
                }
                System.out.println();
                break;
            case 2:
                AlignRange[] results = aligner.align(args[1],args[2]);
                
                for(int i = 0; i < results.length; i++)
                {
                    System.out.println("Hit " + i);
                    System.out.println(args[1].substring(results[i].qStart,results[i].qEnd));
                    System.out.println(args[2].substring(results[i].sStart,results[i].sEnd));
                }
                break;
            case 3:
                int score = 0;
                for(int i = 0; i < args[1].length(); i++)
                {
                    System.out.print(aligner.getScore(args[1].charAt(i),args[2].charAt(i)) + " ");
                    score += aligner.getScore(args[1].charAt(i),args[2].charAt(i));
                }
                System.out.println("Total " + score);
        }
        
    }
    
    public static void usage()
    {
        System.out.println("Usage: java ModuleTest testType query subject");
    }
}