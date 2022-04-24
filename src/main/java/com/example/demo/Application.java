package com.example.demo;

public class Application {

    public static void main(String[] args){
        int size = args.length;
        String helpMessageMain = "###################################\n" +
                             "#   Detect Structural Variants    #\n" +
                             "#       by Kamran Abbaszade       #\n" +
                             "#            24.04.2022           #\n" +
                             "###################################\n\n\n" +
                             "Commands\n"+
                             "-help or -h  \t- printing out help message with instructions\n" +
                             "-start or -s \t- main execution command of application\n";

        String helpMessageStart = "start command\n"+
                                  "Sub-commands\n" +
                "-reference or -r \t- full path to reference FASTA file\n" +
                "-bam or -b       \t- full path to BAM sequence file\n" +
                "Example: -start -reference GRCh38_full_analysis_set_plus_decoy_hla.fa -bam sorted.bam\n";
        if(size == 0 || args[0].contains("-help") || args[0].contains("-h")){
            System.out.println(helpMessageMain);
            System.out.println(helpMessageStart);
        }
        else if(args[0].contains("-start") || args[0].contains("-s")){
            if(((args[1].contains("-reference") || args[1].contains("-r") && args[2].contains(".fa")))){
                if((args[3].contains("-bam") || args[3].contains("-b") && args[4].contains(".bam"))){
                    Algorithm.startWithBam(args[2], args[4]);
                }
                else{
                    System.out.println(helpMessageStart);
                }
            }
        }
        else{
            System.out.println("Unrecognized command.\nRun -help or -h command\n");
        }
    }
}
