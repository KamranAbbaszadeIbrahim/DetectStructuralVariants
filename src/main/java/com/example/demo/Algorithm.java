package com.example.demo;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.*;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Algorithm {
    private static LogUtil log;

    private static class BeginningResult{
        Integer trimLength;
        List<Integer> trimmedCounts;
    }

    private static class GCCorrection{
        List<Double> connectedRD;
        Double globalRD;
    }

    private static class Section{
        List<Integer> startsOfAllSections;
    }

    private static class MergedSection{
        List<Integer> newSectionItems;
    }

    private static class ScanVariants{
        List<Calls> calls;
        List<Boolean> isDeletion;
    }

    private static class Calls{
        int beginning;
        int end;

        Calls(int beginning, int end){
            this.beginning = beginning;
            this.end = end;
        }
    }

    private static class ScanMergedVariants{
        List<Calls> mergedCalls;
        List<Boolean> isDeletionMerged;
    }

    private static class FreezeSectionResults{
        List<Boolean> isFrozen;
        List<Integer> beginningsOfFrozenSections;
    }

    private static class AverageMeanShiftedPartitionResults{
        List<Integer> sectionStart;
        List<Double> newUnfrozenConnectedRD;
    }

    private static class VCF{
        String chromosome;
        Integer position;
        Integer finish;
        String reference;
        String vcfType;
        List<Boolean> isDeletion;

        public VCF(String chromosome, Integer position, Integer finish, String reference, List<Boolean> isDeletion) {
            this.chromosome = chromosome;
            this.position = position;
            this.finish = finish;
            this.reference = reference;
            this.isDeletion = isDeletion;
            if(!isDeletion.isEmpty()){
                this.vcfType = "[DEL]";
            }
            else{
                this.vcfType = "[INS]";
            }
        }

        @Override
        public String toString() {
            return chromosome + " \t" + position.toString() + " \t" + UUID.randomUUID() + " \t" + reference + " \t" + vcfType + " \t.\t_\tSVSIZE=" + (finish - position) + ";";
        }
    }

    private static void init(){
        log = new LogUtil();
    }

    /**
     * startWithBam
     *
     * @param pathToReference = path to reference FASTA file
     * @param pathToBamFile = path to reference BAM file
     */
    public static void startWithBam(String pathToReference, String pathToBamFile) {
        init();
        log.log("LOG_INFO: startWithBam: execution started");
        System.out.println("LOG_INFO: startWithBam: execution started");

        try {
            run(pathToReference, pathToBamFile);
        } catch (Exception e) {
            log.log("LOG_ERROR: startWithBam: " + e.getMessage());
            System.out.println("LOG_ERROR: startWithBam: " + e.getMessage());
        }

        log.log("LOG_INFO: startWithBam: execution completed");
        System.out.println("LOG_INFO: startWithBam: execution completed");
        log.close();
    }

    /**
     * run
     *
     * @param pathToReference = path to reference FASTA file
     * @param pathToBamFile = path to reference BAM file
     */
    private static void run(String pathToReference, String pathToBamFile){
        System.out.println("LOG_INFO: run: execution started");
        log.log("LOG_INFO: run: execution started");

        List<String> chrList = prepareChrList();

        int binLength = 100;

        detectStructuralVariances(pathToBamFile, pathToReference, chrList, binLength);

        log.log("LOG_INFO: run: execution completed");
        System.out.println("LOG_INFO: run: execution completed");
    }

    /**
     * detectStructuralVariances
     *
     * @param pathToBamFile = path to reference BAM file
     * @param pathToReference = path to reference FASTA file
     * @param chrList = ArrayList of chromosome strings
     * @param binLength = size of bin used for calculation of SVs
     */
    private static void detectStructuralVariances(String pathToBamFile, String pathToReference, List<String> chrList, int binLength){
        System.out.println("LOG_INFO: detectStructuralVariances: execution started");
        log.log("LOG_INFO: detectStructuralVariances: execution started");

        SamReader reader = null;

        try{
            reader = SamReaderFactory
                    .makeDefault()
                    .open(new File(pathToBamFile));
        }catch (Exception e){
            log.log("LOG_ERROR: detectStructuralVariances: " + e.getMessage());
            System.out.println("LOG_ERROR: detectStructuralVariances: " + e.getMessage());
        }

        log.log("LOG_INFO: detectStructuralVariances: bamFile opened successfully");
        System.out.println("LOG_INFO: detectStructuralVariances: bamFile opened successfully");

        if(reader == null){
            log.log("LOG_ERROR: detectStructuralVariances: bamFile is null");
            System.out.println("LOG_ERROR: detectStructuralVariances: bamFile is null");
        }

        Map<String, Integer> chromosomeLengths = calculateChromosomeLengths(reader, chrList);

        log.log("LOG_INFO: detectStructuralVariances: chromosomeLengths completed successfully");
        System.out.println("LOG_INFO: detectStructuralVariances: chromosomeLengths completed successfully");

        ReferenceSequenceFile referenceSequenceFile = null;

        try{
            referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(pathToReference), true, true);
        }catch (Exception e){
            log.log("LOG_ERROR: detectStructuralVariances: " + e.getMessage());
            System.out.println("LOG_ERROR: detectStructuralVariances: " + e.getMessage());
        }

        log.log("LOG_INFO: detectStructuralVariances: referenceFile opened successfully");
        System.out.println("LOG_INFO: detectStructuralVariances: referenceFile opened successfully");

        List<VCF> vcfList = new ArrayList<>();

        for(String chr: chrList){
            log.log("LOG_DEBUG: detectStructuralVariances: iteration for " + chr +" started successfully");
            System.out.println("LOG_DEBUG: detectStructuralVariances: iteration for " + chr +" started successfully");

            ReferenceSequence referenceSequence = referenceSequenceFile.getSequence(chr);

            Map<String, List<Integer>> monitoringResults = null;

            try{
                monitoringResults = chromosomeMonitoring(chr, reader, binLength, chromosomeLengths.get(chr), referenceSequence);
                log.log("LOG_INFO: detectStructuralVariances: chromosomeMonitoring: execution completed");
                log.log("LOG_DEBUG: detectStructuralVariances: chromosomeMonitoring: monitoringResults: counts:" + monitoringResults.get("counts"));
                log.log("LOG_DEBUG: detectStructuralVariances: chromosomeMonitoring: monitoringResults: starts:" + monitoringResults.get("starts"));
                log.log("LOG_DEBUG: detectStructuralVariances: chromosomeMonitoring: monitoringResults: gcPercents:" + monitoringResults.get("gcPercents"));
                System.out.println("LOG_INFO: detectStructuralVariances: chromosomeMonitoring: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: chromosomeMonitoring: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: chromosomeMonitoring: " + e.getMessage());
            }

            BeginningResult beginningResult = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: beginning: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: beginning: execution started");
                beginningResult = beginning(monitoringResults.get("counts"));
                log.log("LOG_INFO: detectStructuralVariances: beginning: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: beginning: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: beginning: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: beginning: " + e.getMessage());
            }

            GCCorrection gcCorrection = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: gcCorrection: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: gcCorrection: execution started");
                gcCorrection = gcCorrection(beginningResult.trimmedCounts, monitoringResults.get("gcPercents").subList(beginningResult.trimLength, monitoringResults.get("gcPercents").size()));
                log.log("LOG_INFO: detectStructuralVariances: gcCorrection: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: gcCorrection: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: gcCorrection: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: gcCorrection: " + e.getMessage());
            }

            Section section = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: averageMeanShifting: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: averageMeanShifting: execution started");
                section = averageMeanShifting(beginningResult.trimmedCounts, gcCorrection.connectedRD, gcCorrection.globalRD);
                log.log("LOG_INFO: detectStructuralVariances: averageMeanShifting: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: averageMeanShifting: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: averageMeanShifting: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: averageMeanShifting: " + e.getMessage());
            }

            MergedSection mergedSection = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: sectionMerger: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: sectionMerger: execution started");
                mergedSection = sectionMerger(gcCorrection.connectedRD, section);
                log.log("LOG_INFO: detectStructuralVariances: sectionMerger: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: sectionMerger: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: sectionMerger: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: sectionMerger: " + e.getMessage());
            }

            ScanVariants scanVariants = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: scanVariants: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: scanVariants: execution started");
                scanVariants = scanVariants(mergedSection, gcCorrection.connectedRD);
                log.log("LOG_INFO: detectStructuralVariances: scanVariants: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: scanVariants: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: scanVariants: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: scanVariants: " + e.getMessage());
            }

            ScanMergedVariants scanMergedVariants = null;

            try{
                log.log("LOG_INFO: detectStructuralVariances: scanMergedVariants: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: scanMergedVariants: execution started");
                scanMergedVariants = scanMergedVariants(scanVariants.calls, gcCorrection.connectedRD, scanVariants.isDeletion);
                log.log("LOG_INFO: detectStructuralVariances: scanMergedVariants: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: scanMergedVariants: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: scanMergedVariants: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: scanMergedVariants: " + e.getMessage());
            }

            try{
                log.log("LOG_INFO: detectStructuralVariances: vcfList.addAll:generateVCFCollection: execution started");
                System.out.println("LOG_INFO: detectStructuralVariances: vcfList.addAll:generateVCFCollection: execution started");
                vcfList.addAll(generateVCFCollection(referenceSequenceFile, scanMergedVariants.mergedCalls, binLength, chr, beginningResult.trimLength, scanMergedVariants.isDeletionMerged));
                log.log("LOG_INFO: detectStructuralVariances: vcfList.addAll:generateVCFCollection: execution completed");
                System.out.println("LOG_INFO: detectStructuralVariances: vcfList.addAll:generateVCFCollection: execution completed");
            }catch (Exception e){
                log.log("LOG_ERROR: detectStructuralVariances: vcfList.addAll:generateVCFCollection: " + e.getMessage());
                System.out.println("LOG_ERROR: detectStructuralVariances: vcfList.addAll:generateVCFCollection: " + e.getMessage());
            }
        }

        try{
            log.log("LOG_INFO: detectStructuralVariances: writeVCFFile: execution started");
            System.out.println("LOG_INFO: detectStructuralVariances: writeVCFFile: execution started");
            writeVCFFile("DetectedStructuralVariants.vcf", vcfList);
            log.log("LOG_INFO: detectStructuralVariances: writeVCFFile: execution completed");
            System.out.println("LOG_INFO: detectStructuralVariances: writeVCFFile: execution completed");
        }catch (Exception e){
            log.log("LOG_ERROR: detectStructuralVariances: writeVCFFile: " + e.getMessage());
            System.out.println("LOG_ERROR: detectStructuralVariances: writeVCFFile: " + e.getMessage());
        }

        log.log("LOG_INFO: detectStructuralVariances: execution completed");
        System.out.println("LOG_INFO: detectStructuralVariances: execution completed");
    }

    /**
     * calculateChromosomeLengths
     *
     * @param reader = SamReader file containing link to BAM file
     * @param chrList = ArrayList of chromosome strings
     * @return Map<String, Integer>
     */
    private static Map<String, Integer> calculateChromosomeLengths(SamReader reader, List<String> chrList){
        System.out.println("LOG_INFO: calculateChromosomeLengths: execution started");
        log.log("LOG_INFO: calculateChromosomeLengths: execution started");

        List<SAMSequenceRecord> sequenceRecords = reader.getFileHeader().getSequenceDictionary().getSequences();
        Map<String, Integer> chromosomeLengths = new HashMap<>();

        for(SAMSequenceRecord record: sequenceRecords){
            for(String chr: chrList){
                if(chr.equals(record.getSequenceName())){
                    chromosomeLengths.put(record.getSequenceName(), record.getSequenceLength());
                }
            }
        }

        if(chromosomeLengths.isEmpty()){
            log.log("LOG_ERROR: calculateChromosomeLengths: chromosomeLengths is empty");
            System.out.println("LOG_ERROR: calculateChromosomeLengths: chromosomeLengths is empty");
        }

        log.log("LOG_DEBUG: calculateChromosomeLengths: chromosomeLengths: " + chromosomeLengths);
        System.out.println("LOG_DEBUG: calculateChromosomeLengths: chromosomeLengths: " + chromosomeLengths);

        log.log("LOG_INFO: calculateChromosomeLengths: execution completed");
        System.out.println("LOG_INFO: calculateChromosomeLengths: execution completed");
        return chromosomeLengths;
    }

    /**
     * chromosomeMonitoring
     *
     * @param chr = chromosome string
     * @param reader = SamReader object containing link to BAM file
     * @param binLength = length of bin
     * @param chromosomeLength = length of current chromosome
     * @param referenceSequence = reference sequence for current chromosome from FASTA file
     * @return Map<String, List<Integer>>
     */
    private static Map<String, List<Integer>> chromosomeMonitoring(String chr, SamReader reader, int binLength, int chromosomeLength, ReferenceSequence referenceSequence){
        System.out.println("LOG_INFO: chromosomeMonitoring: execution started");
        System.out.println("LOG_DEBUG: chromosomeMonitoring: chr="+chr+"; binLength="+binLength+"; chrLength="+chromosomeLength);
        log.log("LOG_INFO: chromosomeMonitoring: execution started");
        log.log("LOG_DEBUG: chromosomeMonitoring: chr="+chr+"; binLength="+binLength+"; chrLength="+chromosomeLength);

        Map<String, List<Integer>> monitoringResult = new HashMap<>();

        List<Integer> starts = new ArrayList<>();
        List<Integer> counts = new ArrayList<>();
        List<Integer> gcPercents = new ArrayList<>();

        int start = 0;
        int partition = 0;

        while(start < chromosomeLength){
            if(chromosomeLength - start < binLength){
                break;
            }

            SAMRecordIterator reads = reader.queryContained(chr, start, start + binLength);
            int i = reads.toList().size();

            String fragment = referenceSequence.getBaseString().substring(start, start + binLength);


            double gcPercent = calculateGCContent(fragment);

            start = start + binLength;

            counts.add(i);
            starts.add(start);
            gcPercents.add((int)gcPercent);

            partition = partition + 1;

            if(partition == 1000){
                double progress = ((double) start * 100)/chromosomeLength;
                System.out.println(String.format("Progress for %s --> %.3f%s", chr, progress, "%"));
                log.log(String.format("Progress for %s --> %.3f%s", chr, progress, "%"));

                partition = 0;
            }

        }

        monitoringResult.put("counts", counts);
        monitoringResult.put("starts", starts);
        monitoringResult.put("gcPercents", gcPercents);

        System.out.println("LOG_INFO: chromosomeMonitoring: execution completed");
        log.log("LOG_INFO: chromosomeMonitoring: execution completed");

        return monitoringResult;
    }

    /**
     * beginning
     *
     * @param counts = scanned chromosome counts
     * @return BeginningResult
     */
    private static BeginningResult beginning(List<Integer> counts){

        BeginningResult beginningResult = new BeginningResult();
        int i = 0;

        for(Integer read: counts){
            if(read == 0){
                i++;
            }
            else{
                break;
            }
        }

        beginningResult.trimLength = i;
        beginningResult.trimmedCounts = counts.subList(i, counts.size());

        return beginningResult;
    }

    /**
     * gcCorrection
     *
     * @param trimmedCounts
     * @param gcPercents
     * @return GCCorrection
     */
    private static GCCorrection gcCorrection(List<Integer> trimmedCounts, List<Integer> gcPercents){
        GCCorrection gcCorrection = new GCCorrection();
        List<Integer> uniqueGC = gcPercents.stream().distinct().collect(Collectors.toList());
        Map<Integer, Double> readGCS = new HashMap<>();

        for(Integer element: uniqueGC){
            List<Integer> read = new ArrayList<>();
            for(int i=0; i<gcPercents.size(); i++){
                if(Objects.equals(gcPercents.get(i), element)){
                    read.add(trimmedCounts.get(i));
                }
            }
            double rdLocalMean = 0;

            for(int i=0; i<read.size(); i++){
                rdLocalMean = rdLocalMean + read.get(i);
            }
            rdLocalMean = rdLocalMean / read.size();

            if(rdLocalMean != 0){
                readGCS.put(element, rdLocalMean);
            }
            else{
                readGCS.put(element, 1.0);
            }
        }

        double rdGlobalMean = 0;
        for(int i=0; i<trimmedCounts.size(); i++){
            rdGlobalMean = rdGlobalMean + trimmedCounts.get(i);
        }
        rdGlobalMean = rdGlobalMean / trimmedCounts.size();

        List<Double> readGCSMix = new ArrayList<>();
        for(Integer percent: gcPercents){
            readGCSMix.add(readGCS.get(percent));
        }

        List<Double> readAdjusted = new ArrayList<>();
        for(int i = 0; i<gcPercents.size(); i++){
            readAdjusted.add((rdGlobalMean / readGCSMix.get(i)) * trimmedCounts.get(i));
        }

        gcCorrection.connectedRD = readAdjusted;
        gcCorrection.globalRD = rdGlobalMean;

        return gcCorrection;
    }

    /**
     * calculateGCContent
     *
     * @param fragment = fragment of sequence
     * @return Double
     */
    private static Double calculateGCContent(String fragment){

        double gc = 0;
        for(Character c: fragment.toCharArray()){
            switch (c){
                case 'G':
                case 'C':
                case 'g':
                case 'c':
                case 'S':
                case 's': gc++;
                break;
            }
        }


        return gc * 100 / fragment.length();
    }

    /**
     * averageMeanShifting
     *
     * @param trimmedCounts
     * @param connectedRD
     * @param globalRD
     * @return Section
     */
    private static Section averageMeanShifting(List<Integer> trimmedCounts, List<Double> connectedRD, Double globalRD){
        System.out.println("LOG_INFO: averageMeanShifting: execution started");
        log.log("LOG_INFO: averageMeanShifting: execution started");

        Section section = null;
        double mean = MathUtility.meanDouble(connectedRD);
        double standardDeviation = MathUtility.standardDeviationDouble(connectedRD);

        List<Double> hyperParam = new ArrayList<>();
        for(int i=0; i< connectedRD.size(); i++){
            if(connectedRD.get(i) > globalRD / 4){
                hyperParam.add(Math.sqrt(trimmedCounts.get(i) / globalRD) * standardDeviation);
            }
            else{
                hyperParam.add(standardDeviation/2);
            }
        }

        int hyperBatch = 2;
        int limit = 128;

        List<Boolean> isFrozen = new ArrayList<>(Arrays.asList(new Boolean[trimmedCounts.size()]));
        Collections.fill(isFrozen, Boolean.FALSE);

        List<Integer> beginningsOfFrozenSections = new ArrayList<>();

        while(hyperBatch < limit) {
            section = averageMeanShiftingStep(connectedRD, hyperBatch, hyperParam, isFrozen, beginningsOfFrozenSections);
            FreezeSectionResults freezeSectionResults = freezeSections(section, connectedRD, hyperParam);
            isFrozen = freezeSectionResults.isFrozen;
            beginningsOfFrozenSections = freezeSectionResults.beginningsOfFrozenSections;
            hyperBatch = hyperBatch + 1;
        }

        System.out.println("LOG_INFO: averageMeanShifting: execution completed");
        log.log("LOG_INFO: averageMeanShifting: execution completed");
        return section;
    }

    /**
     * averageMeanShiftingStep
     *
     * @param connectedRD
     * @param hyperBatch
     * @param hyperParam
     * @param isFrozen
     * @param beginningsOfFrozenSections
     * @return section
     */
    private static Section averageMeanShiftingStep(List<Double> connectedRD, Integer hyperBatch, List<Double> hyperParam, List<Boolean> isFrozen, List<Integer> beginningsOfFrozenSections){
        System.out.println("LOG_INFO: averageMeanShiftingStep: execution started");
        log.log("LOG_INFO: averageMeanShiftingStep: execution started");

        List<Double> unfrozenConnectedRD = new ArrayList<>();
        List<Double> unfrozenHyperParam = new ArrayList<>();

        for(int i=0; i< connectedRD.size(); i++){
            if(!isFrozen.get(i)){
                unfrozenConnectedRD.add(connectedRD.get(i));
            }
        }
        for(int i=0; i< hyperParam.size(); i++){
            if(!isFrozen.get(i)){
                unfrozenHyperParam.add(hyperParam.get(i));
            }
        }

        List<Double> shiftedRD = new ArrayList<>(unfrozenConnectedRD);
        AverageMeanShiftedPartitionResults partitionResults = new AverageMeanShiftedPartitionResults();
        partitionResults.newUnfrozenConnectedRD = shiftedRD;
        partitionResults.sectionStart = new ArrayList<>();

        for(int i=0; i< 3; i++){
            List<Double> gradientResult = averageMeanShiftingGradient(partitionResults.newUnfrozenConnectedRD, hyperBatch, unfrozenHyperParam);
            partitionResults = averageMeanShiftingPartition(partitionResults.newUnfrozenConnectedRD, gradientResult);
        }

        System.out.println("LOG_INFO: averageMeanShiftingStep: execution completed");
        log.log("LOG_INFO: averageMeanShiftingStep: execution completed");

        return frozenSectionInjection(partitionResults.sectionStart, beginningsOfFrozenSections, isFrozen);
    }

    /**
     * frozenSectionInjection
     *
     * @param sectionStart
     * @param beginningsOfFrozenSections
     * @param isFrozen
     * @return section
     */
    private static Section frozenSectionInjection(List<Integer> sectionStart, List<Integer> beginningsOfFrozenSections, List<Boolean> isFrozen){
        System.out.println("LOG_INFO: frozenSectionInjection: execution started");
        log.log("LOG_INFO: frozenSectionInjection: execution started");

        List<Integer> startsOfAllSections = new ArrayList<>();
        int shift = 0;
        int indexOfUnfrozen = 0;
        int indexOfFrozen = 0;

        for(int i = 0; i< isFrozen.size(); i++){
            if(isFrozen.get(i)){
                shift = shift + 1;

                if(indexOfFrozen < beginningsOfFrozenSections.size() && beginningsOfFrozenSections.get(indexOfFrozen) == i){
                    startsOfAllSections.add(i);
                    indexOfFrozen = indexOfUnfrozen + 1;
                }
            }
            else{
                if(indexOfUnfrozen < sectionStart.size() && sectionStart.get(indexOfUnfrozen) + shift == i){
                    startsOfAllSections.add(i);
                    indexOfUnfrozen = indexOfUnfrozen + 1;
                }
            }
        }

        Section section = new Section();
        section.startsOfAllSections = startsOfAllSections;
        System.out.println("LOG_INFO: frozenSectionInjection: execution completed");
        log.log("LOG_INFO: frozenSectionInjection: execution completed");
        return section;
    }

    /**
     * freezeSections
     *
     * @param section
     * @param connectedRD
     * @param hyperParam
     * @return freezeSectionResults
     */
    private static FreezeSectionResults freezeSections(Section section, List<Double> connectedRD, List<Double> hyperParam){
        System.out.println("LOG_INFO: freezeSections: execution started");
        log.log("LOG_INFO: freezeSections: execution started");

        double meanRD = MathUtility.meanDouble(connectedRD);

        List<Boolean> isFrozen = new ArrayList<>(Arrays.asList(new Boolean[connectedRD.size()]));
        Collections.fill(isFrozen, Boolean.FALSE);

        List<Integer> beginningsOfFrozenSections = new ArrayList<>();

        int n = section.startsOfAllSections.size();

        List<Double> means = new ArrayList<>(Arrays.asList(new Double[n]));
        Collections.fill(means, Double.MIN_VALUE);

        List<Double> sides = new ArrayList<>(Arrays.asList(new Double[n]));
        Collections.fill(sides, Double.MIN_VALUE);

        List<Double> names = new ArrayList<>(Arrays.asList(new Double[n]));
        Collections.fill(names, Double.MIN_VALUE);

        List<Integer> sectionEdges = section.startsOfAllSections;

        sectionEdges.add(n);

        int lengthOfGenome = connectedRD.size();

        for(int i = 0; i< n; i++){
            means.set(i, MathUtility.meanDouble(connectedRD.subList(sectionEdges.get(i), sectionEdges.get(i+1))));
            sides.set(i, MathUtility.standardDeviationDouble(connectedRD.subList(sectionEdges.get(i), sectionEdges.get(i+1))));
            names.set(i, (double) (sectionEdges.get(i + 1) - sectionEdges.get(i)));
        }

        for(int i = 0; i< n; i++){
            double toggle = ((meanRD - means.get(i)) / sides.get(i)) * Math.sqrt(names.get(i));

            double particle = 2 * (MathUtility.studentTCDF(-Math.abs(toggle), names.get(i) - 1));

            int sectionLength = sectionEdges.get(i + 1) - sectionEdges.get(i);

            double particleCoordinate = (particle * 0.99 * lengthOfGenome) / sectionLength;

            boolean checkPoint1 = particleCoordinate < 0.05;
            boolean checkPoint2;
            boolean checkPoint3;

            if(i > 0){
                double meanVariance = means.get(i - 1) - means.get(i);
                double side1 = sides.get(i - 1);
                double side2 = sides.get(i);
                double node1 = sectionEdges.get(i) - sectionEdges.get(i - 1);
                double node2 = sectionLength;
                double tumbler = meanVariance / Math.sqrt(Math.pow(side1,2) / node1 + Math.pow(side2, 2) / node2);
                double newParticle = 2 * (MathUtility.studentTCDF(-Math.abs(tumbler), 1));
                double newParticleCoordinate = (newParticle * 0.99 * lengthOfGenome) / (node1 + node2);
                checkPoint2 = newParticleCoordinate < 0.01 || Math.abs(meanVariance) >= 2 * hyperParam.get(i);
            }
            else{
                checkPoint2 = true;
            }
            if(i < n - 1){
                double meanVariance = means.get(i + 1) - means.get(i);
                double side1 = sides.get(i + 1);
                double side2 = sides.get(i);
                double node1 = sectionEdges.get(i + 1) - sectionEdges.get(i);
                double node2 = sectionLength;
                double tag = meanVariance / Math.sqrt(Math.pow(side1, 2) / node1 + Math.pow(side2, 2) / node2);
                double anotherParticle = 2 * (MathUtility.studentTCDF(-Math.abs(tag), 1));
                double anotherParticleCoordinate = (anotherParticle * 0.99 * lengthOfGenome) / (node1 + node2);
                checkPoint3 = anotherParticleCoordinate < 0.01 || Math.abs(meanVariance) >= 2 * hyperParam.get(i);
            }
            else{
                checkPoint3 = true;
            }

            if(checkPoint1 && checkPoint2 && checkPoint3){
                beginningsOfFrozenSections.add(section.startsOfAllSections.get(i));
                for(int x = sectionEdges.get(i); x < sectionEdges.get(i + 1); x++){
                    isFrozen.set(x, Boolean.TRUE);
                }
            }
        }

        FreezeSectionResults freezeSectionResults = new FreezeSectionResults();
        freezeSectionResults.beginningsOfFrozenSections = beginningsOfFrozenSections;
        freezeSectionResults.isFrozen = isFrozen;
        System.out.println("LOG_INFO: freezeSections: execution completed");
        log.log("LOG_INFO: freezeSections: execution completed");
        return freezeSectionResults;
    }

    /**
     * averageMeanShiftingGradient
     *
     * @param connectedRD
     * @param hyperBatch
     * @param hyperParam
     * @return List<Double>
     */
    private static List<Double> averageMeanShiftingGradient(List<Double> connectedRD, Integer hyperBatch, List<Double> hyperParam){
        System.out.println("LOG_INFO: averageMeanShiftingGradient: execution started");
        log.log("LOG_INFO: averageMeanShiftingGradient: execution started");

        Integer start = 0;
        Integer maximalLength = 40_000_000;
        List<Double> gradient = new ArrayList<>(Arrays.asList(new Double[connectedRD.size()]));
        Collections.fill(gradient, Double.MIN_VALUE);

        while(start < connectedRD.size()){
            int n = Math.min(connectedRD.size()-start, maximalLength);
            List<List<List<Integer>>> indexes = buildIndexMatrix(n, hyperBatch);
            indexes = adjustIndexMatrix(indexes, start, "+");

            List<List<Integer>> rightAdjustedIndexMatrix = adjustRightAngle(indexes);
            correctRightAdjustedAngle(indexes, rightAdjustedIndexMatrix, connectedRD);

            List<List<Integer>> leftAdjustedIndexMatrix = adjustLeftAngle(indexes);
            correctLeftAdjustedAngle(indexes, leftAdjustedIndexMatrix, connectedRD);

            List<List<Integer>> mergedIndexMatrix = mergeMatrices(indexes.get(0), indexes.get(0));
            List<List<Integer>> mergedLeftRightAngles = mergeMatrices(leftAdjustedIndexMatrix, rightAdjustedIndexMatrix);
            List<List<Integer>> subtractedAngles = subtractionOfMatrices(mergedIndexMatrix, mergedLeftRightAngles);

            List<List<Double>> exponent1 = exponent1(subtractedAngles, hyperBatch);
            List<List<Double>> exponent2 = exponent2(mergedIndexMatrix, mergedLeftRightAngles, connectedRD, hyperParam);

            gradient = updateGradient(gradient, multiplier(subtractedAngles, exponent1, exponent2), start, n);
            start = n;
        }

        System.out.println("LOG_INFO: averageMeanShiftingGradient: execution completed");
        log.log("LOG_INFO: averageMeanShiftingGradient: execution completed");
        return gradient;
    }

    private static List<List<List<Integer>>> buildIndexMatrix(int n, int m){
        List<List<List<Integer>>> indexes = new ArrayList<>();

        {
            List<List<Integer>> rec = new ArrayList<>();
            for(int j=0; j<n; j++){
                List<Integer> elem = new ArrayList<>(Arrays.asList(new Integer[m]));
                Collections.fill(elem, j);
                rec.add(elem);
            }
            indexes.add(rec);
        }
        {
            List<Integer> elem = new ArrayList<>();
            for(int j=0; j<m; j++){
                elem.add(j);
            }
            List<List<Integer>> rec = new ArrayList<>();
            for(int j=0; j< n; j++){
                rec.add(elem);
            }
            indexes.add(rec);
        }
        return indexes;
    }

    private static List<List<List<Integer>>> adjustIndexMatrix(List<List<List<Integer>>> indexMatrix, int adjustment, String sign){
        if(sign.contains("+")){
            for(int k=0; k < indexMatrix.size(); k++){
                for(int i=0; i< indexMatrix.get(k).size(); i++){
                    List<Integer> tmp = new ArrayList<>();
                    for(int j=0; j < indexMatrix.get(k).get(i).size(); j++){
                        tmp.add(indexMatrix.get(k).get(i).get(j) + adjustment);
                    }
                    indexMatrix.get(k).set(i, tmp);
                }
            }
        }
        else if(sign.contains("-")){
            for(int k=0; k < indexMatrix.size(); k++){
                for(int i=0; i< indexMatrix.get(k).size(); i++){
                    List<Integer> tmp = new ArrayList<>();
                    for(int j=0; j < indexMatrix.get(k).get(i).size(); j++){
                        tmp.add(indexMatrix.get(k).get(i).get(j) - adjustment);
                    }
                    indexMatrix.get(k).set(i, tmp);
                }
            }
        }
        else{
            System.out.println("Can not apply provided operation on index matrix");
        }
        return indexMatrix;
    }

    private static List<List<Integer>> adjustRightAngle(List<List<List<Integer>>> indexMatrix){
        List<List<Integer>> rightAdjustedMatrix = new ArrayList<>();
        for(int i = 0; i< indexMatrix.get(0).size(); i++){
            List<Integer> tmp = new ArrayList<>();
            for(int j = 0; j< indexMatrix.get(0).get(i).size(); j++){
                tmp.add(indexMatrix.get(1).get(i).get(j) + indexMatrix.get(0).get(i).get(j) + 1);
            }
            rightAdjustedMatrix.add(tmp);
        }

        return rightAdjustedMatrix;
    }

    private static void correctRightAdjustedAngle(List<List<List<Integer>>> indexMatrix, List<List<Integer>> rightAdjustedMatrix, List<Double> connectedRD){
        for(int i=0; i<rightAdjustedMatrix.size(); i++){
            for(int j=0; j<rightAdjustedMatrix.get(i).size(); j++){
                if(rightAdjustedMatrix.get(i).get(j) >= connectedRD.size()){
                    List<Integer> tmp = rightAdjustedMatrix.get(i);
                    tmp.set(j, indexMatrix.get(0).get(i).get(j));
                    rightAdjustedMatrix.set(i, tmp);
                }
            }
        }
    }

    private static List<List<Integer>> adjustLeftAngle(List<List<List<Integer>>> indexMatrix){
        List<List<Integer>> rightAdjustedMatrix = new ArrayList<>();
        for(int i = 0; i< indexMatrix.get(0).size(); i++){
            List<Integer> tmp = new ArrayList<>();
            for(int j = 0; j< indexMatrix.get(0).get(i).size(); j++){
                tmp.add(indexMatrix.get(0).get(i).get(j) - indexMatrix.get(1).get(i).get(j) - 1);
            }
            rightAdjustedMatrix.add(tmp);
        }

        return rightAdjustedMatrix;
    }

    private static void correctLeftAdjustedAngle(List<List<List<Integer>>> indexMatrix, List<List<Integer>> leftAdjustedMatrix, List<Double> connectedRD){
        for(int i=0; i<leftAdjustedMatrix.size(); i++){
            for(int j=0; j<leftAdjustedMatrix.get(i).size(); j++){
                if(leftAdjustedMatrix.get(i).get(j) < 0){
                    List<Integer> tmp = leftAdjustedMatrix.get(i);
                    tmp.set(j, indexMatrix.get(0).get(i).get(j));
                    leftAdjustedMatrix.set(i, tmp);
                }
            }
        }
    }

    private static List<List<Integer>> mergeMatrices(List<List<Integer>> a, List<List<Integer>> b){
        List<List<Integer>> result = new ArrayList<>();
        for(int i=0; i<a.size(); i++){
            List<Integer> tmp = new ArrayList<>();
            tmp.addAll(a.get(i));
            tmp.addAll(b.get(i));
            result.add(tmp);
        }
        return result;
    }

    private static List<List<Integer>> subtractionOfMatrices(List<List<Integer>> a, List<List<Integer>> b){
        List<List<Integer>> result = new ArrayList<>();
        for(int i=0; i<a.size(); i++){
            List<Integer> tmp = new ArrayList<>();
            for(int j=0; j<a.get(i).size(); j++){
                tmp.add(b.get(i).get(j) - a.get(i).get(j));
            }
            result.add(tmp);
        }
        return result;
    }

    private static List<List<Double>> exponent1(List<List<Integer>> a, Integer hyperBatch){
        List<List<Double>> result = new ArrayList<>();
        for(int i=0; i<a.size(); i++){
            List<Double> tmp = new ArrayList<>();
            for(int j=0; j<a.get(i).size(); j++){
                double square = Math.pow(a.get(i).get(j), 2);
                double inverse = -square;
                double divider = 2 * Math.pow(hyperBatch, 2);
                double value = Math.exp(inverse/divider);
                tmp.add(value);
            }
            result.add(tmp);
        }
        return result;
    }

    private static List<List<Double>> exponent2(List<List<Integer>> mergedIndexMatrix, List<List<Integer>> mergedLeftRightAngles, List<Double> connectedRD, List<Double> hyperParam){
        List<List<Double>> result = new ArrayList<>();

        for(int i = 0; i<mergedIndexMatrix.size(); i++){
            List<Double> tmp = new ArrayList<>();
            for(int j = 0; j<mergedIndexMatrix.get(i).size(); j++){
                double x = connectedRD.get(mergedLeftRightAngles.get(i).get(j));
                double y = connectedRD.get(mergedIndexMatrix.get(i).get(j));
                double z = x - y;
                double square = Math.pow(z, 2);
                double inverse = -square;
                double newHyperParam = hyperParam.get(mergedIndexMatrix.get(i).get(j));
                double squareHyperParam = Math.pow(newHyperParam, 2);
                double divider = 2 * squareHyperParam;
                double value = Math.exp(inverse/divider);
                tmp.add(value);
            }
            result.add(tmp);
        }

        return result;
    }

    private static List<List<Double>> multiplier(List<List<Integer>> subtractedAngles, List<List<Double>> exponent1, List<List<Double>> exponent2){
        List<List<Double>> result = new ArrayList<>();
        for(int i = 0; i< subtractedAngles.size(); i++){
            List<Double> tmp = new ArrayList<>();
            for(int j=0; j< subtractedAngles.get(i).size(); j++){
                double x = subtractedAngles.get(i).get(j) * exponent1.get(i).get(j) * exponent2.get(i).get(j);
                tmp.add(x);
            }
            result.add(tmp);
        }

        return result;
    }

    private static List<Double> updateGradient(List<Double> gradient, List<List<Double>> multiplier, Integer start, Integer n){
        List<Double> sums = new ArrayList<>();

        for(int i=0; i< multiplier.size(); i++){
            double sum = 0;
            for(int j=0; j<multiplier.get(i).size(); j++){
                sum = sum + multiplier.get(i).get(j);
            }
            sums.add(sum);
        }

        for(int i = start, x = 0; i< (n + start); i++, x++){
            gradient.set(start, sums.get(x));
        }

        return gradient;
    }

    /**
     * averageMeanShiftingPartition
     *
     * @param unfrozenConnectedRD
     * @param gradientResult
     * @return AverageMeanShiftedPartitionResults
     */
    private static AverageMeanShiftedPartitionResults averageMeanShiftingPartition(List<Double> unfrozenConnectedRD, List<Double> gradientResult){
        System.out.println("LOG_INFO: averageMeanShiftingPartition: execution started");
        log.log("LOG_INFO: averageMeanShiftingPartition: execution started");
        AverageMeanShiftedPartitionResults partitionResults = new AverageMeanShiftedPartitionResults();

        Integer start = 0;
        List<Integer> sectionStart = new ArrayList<>();
        sectionStart.add(0);
        List<Double> newUnfrozenConnectedRD = new ArrayList<>(Arrays.asList(new Double[unfrozenConnectedRD.size()]));
        Collections.fill(newUnfrozenConnectedRD, Double.MIN_VALUE);


        for(int i=0; i < gradientResult.size(); i++){
            if(i == gradientResult.size()-1){
                double mean = MathUtility.meanDouble(unfrozenConnectedRD);
                if((start + (i+1)) < unfrozenConnectedRD.size()){
                    for(int j = start; j < (start + (i+1)); j++){
                        newUnfrozenConnectedRD.set(j, mean);
                    }
                }
            }
            else if(gradientResult.get(i) <= 0 && gradientResult.get(i) < gradientResult.get(i+1) && 0 < gradientResult.get(i+1)){
                double mean = MathUtility.meanDouble(unfrozenConnectedRD);;

                if((start + (i+1)) < unfrozenConnectedRD.size()){
                    for(int j = start; j < (start + (i+1)); j++){
                        newUnfrozenConnectedRD.set(j, mean);
                    }
                }

                start = i + 1;
                sectionStart.add(start);
            }
        }

        partitionResults.sectionStart = sectionStart;
        partitionResults.newUnfrozenConnectedRD = newUnfrozenConnectedRD;
        System.out.println("LOG_INFO: averageMeanShiftingPartition: execution completed");
        log.log("LOG_INFO: averageMeanShiftingPartition: execution completed");
        return partitionResults;
    }

    /**
     * prepareChrList
     *
     * @return List<String>
     */
    private static List<String> prepareChrList(){
        System.out.println("LOG_INFO: prepareChrList: execution started");
        log.log("LOG_INFO: prepareChrList: execution started");

        List<String> chrList = new ArrayList<>();
        for(int i = 1; i<=22; i++){
            chrList.add("chr"+i);
            if(i==22){
                chrList.add("chrX");
            }
        }

        System.out.println("LOG_INFO: prepareChrList: execution completed");
        log.log("LOG_INFO: prepareChrList: execution completed");
        return chrList;
    }

    /**
     * sectionMerger
     *
     * @param connectedRD
     * @param section
     * @return MergedSection
     */
    private static MergedSection sectionMerger(List<Double> connectedRD, Section section){
        System.out.println("LOG_INFO: sectionMerger: execution started");
        log.log("LOG_INFO: sectionMerger: execution started");
        double threshold = MathUtility.meanDouble(connectedRD) / 4;
        List<Integer> newSectionItems = new ArrayList<>();
        newSectionItems.add(0);
        section.startsOfAllSections.add(connectedRD.size());
        for(int i = 2; i< section.startsOfAllSections.size(); i++){
            double lastMean = MathUtility.meanDouble(connectedRD.subList(newSectionItems.get(newSectionItems.size()-1), section.startsOfAllSections.get(i-1)));

            double currentMean = MathUtility.meanDouble(connectedRD.subList(section.startsOfAllSections.get(i-1), section.startsOfAllSections.get(i)));

            if(Math.abs(lastMean - currentMean) > threshold){
                newSectionItems.add(section.startsOfAllSections.get(i - 1));
            }
        }

        MergedSection mergedSection = new MergedSection();
        mergedSection.newSectionItems = newSectionItems;
        System.out.println("LOG_INFO: sectionMerger: execution completed");
        log.log("LOG_INFO: sectionMerger: execution completed");
        return mergedSection;
    }

    /**
     * scanVariants
     *
     * @param section
     * @param connectedRD
     * @return ScanVariants
     */
    private static ScanVariants scanVariants(MergedSection section, List<Double> connectedRD){
        System.out.println("LOG_INFO: scanVariants: execution started");
        log.log("LOG_INFO: scanVariants: execution started");
        section.newSectionItems.add(connectedRD.size());

        double globalRD = MathUtility.meanDouble(connectedRD);

        List<Calls> calls = new ArrayList<>();

        List<Boolean> isDeletion = new ArrayList<>();

        for(int i=1; i< section.newSectionItems.size(); i++){
            int beginning = section.newSectionItems.get(i - 1);

            int end = section.newSectionItems.get(i);

            double radialSection = MathUtility.meanDouble(connectedRD.subList(beginning, end));

            double stockSection = MathUtility.standardDeviationDouble(connectedRD.subList(beginning, end));

            int n = end - beginning;

            if(Math.abs(globalRD - radialSection) < globalRD / 4){
                continue;
            }

            if(stockSection == 0){
                calls.add(new Calls(beginning, end));
                isDeletion.add((globalRD - radialSection) > 0);
                continue;
            }

            double tStatistics = ((globalRD - radialSection) / stockSection) * Math.sqrt(n);

            double particle = 2 * (MathUtility.studentTCDF(-Math.abs(tStatistics), n - 1));

            double particleCoordinate = (particle * 0.99 * connectedRD.size()) / n;

            if(particleCoordinate < 0.0005){
                calls.add(new Calls(beginning, end));
                isDeletion.add((globalRD - radialSection) > 0);
            }
        }

        ScanVariants scanVariants = new ScanVariants();
        scanVariants.calls = calls;
        scanVariants.isDeletion = isDeletion;

        System.out.println("LOG_INFO: scanVariants: execution completed");
        log.log("LOG_INFO: scanVariants: execution completed");
        return scanVariants;
    }


    /**
     * scanMergedVariants
     *
     * @param calls
     * @param connectedRD
     * @param isDeletion
     * @return ScanMergedVariants
     */
    private static ScanMergedVariants scanMergedVariants(List<Calls> calls, List<Double> connectedRD, List<Boolean> isDeletion){
        System.out.println("LOG_INFO: scanMergedVariants: execution started");
        log.log("LOG_INFO: scanMergedVariants: execution started");
        List<Calls> mergedCalls = new ArrayList<>();

        List<Boolean> isDeletionMerged = new ArrayList<>();

        int lengthOfGenome = connectedRD.size();

        for(int i=1; i< calls.size(); i++){
            List<Double> firstScannedCall = connectedRD.subList(calls.get(i - 1).beginning, calls.get(i - 1).end);

            List<Double> secondScannedCall = connectedRD.subList(calls.get(i).beginning, calls.get(i).end);

            List<Double> region = connectedRD.subList(calls.get(i-1).end, calls.get(i).beginning);

            double firstMean = MathUtility.meanDouble(firstScannedCall);

            double secondMean = MathUtility.meanDouble(secondScannedCall);

            double regionMean = MathUtility.meanDouble(region);

            double firstStandardDeviation = MathUtility.standardDeviationDouble(firstScannedCall);

            double secondStandardDeviation = MathUtility.standardDeviationDouble(secondScannedCall);

            double regionStandardDeviation = MathUtility.standardDeviationDouble(region);

            int n1 = firstScannedCall.size();

            int n2 = secondScannedCall.size();

            int nr = region.size();

            double firstToggle = (firstMean - regionMean) / Math.sqrt(Math.pow(firstStandardDeviation, 2) / n1 + Math.pow(regionStandardDeviation, 2) / nr);

            double secondToggle = (secondMean - regionMean) / Math.sqrt(Math.pow(secondStandardDeviation, 2) / n2 + Math.pow(regionStandardDeviation, 2) / nr);

            double firstParticle = 2 * (MathUtility.studentTCDF(-Math.abs(firstToggle), 1));

            double secondParticle = 2 * (MathUtility.studentTCDF(-Math.abs(secondToggle), 1));

            double firstParticleCoordinate = (firstParticle * 0.01 * lengthOfGenome) / (n1 + nr);

            double secondParticleCoordinate = (secondParticle * 0.01 * lengthOfGenome) / (n2 + nr);

            if(firstParticleCoordinate > 0.01 && secondParticleCoordinate > 0.01){
                calls.set(i, new Calls(calls.get(i - 1).beginning, calls.get(i).beginning));
            }
            else{
                mergedCalls.add(new Calls(calls.get(i - 1).beginning, calls.get(i - 1).end));
                isDeletionMerged.add(isDeletion.get(i - 1));
            }
        }

        mergedCalls.add(calls.get(calls.size() - 1));
        isDeletionMerged.add(isDeletion.get(isDeletion.size() - 1));

        ScanMergedVariants scanMergedVariants = new ScanMergedVariants();
        scanMergedVariants.mergedCalls = mergedCalls;
        scanMergedVariants.isDeletionMerged = isDeletionMerged;

        System.out.println("LOG_INFO: scanMergedVariants: execution completed");
        log.log("LOG_INFO: scanMergedVariants: execution completed");
        return scanMergedVariants;
    }

    /**
     * generateVCFCollection
     *
     * @param referenceSequenceFile
     * @param calls
     * @param binLength
     * @param chromosome
     * @param trimLength
     * @param isDeletion
     * @return List<VCF>
     */
    private static List<VCF> generateVCFCollection(ReferenceSequenceFile referenceSequenceFile, List<Calls> calls, int binLength, String chromosome, Integer trimLength, List<Boolean> isDeletion){
        System.out.println("LOG_INFO: generateVCFCollection: execution started");
        log.log("LOG_INFO: generateVCFCollection: execution started");
        List<VCF> vcfList = new ArrayList<>();
        ReferenceSequence referenceSequence = referenceSequenceFile.getSequence(chromosome);
        for(Calls call: calls){
            int start = call.beginning * binLength + trimLength;
            int end = call.end * binLength + trimLength;
            String ref = String.valueOf((char)referenceSequence.getBases()[start]);
            VCF vcf = new VCF(chromosome,start, end, ref, isDeletion);
            vcfList.add(vcf);
        }
        System.out.println("LOG_INFO: generateVCFCollection: execution completed");
        log.log("LOG_INFO: generateVCFCollection: execution completed");
        return vcfList;
    }

    /**
     * writeVCFFile
     *
     * @param pathToVCFFile
     * @param vcfList
     */
    private static void writeVCFFile(String pathToVCFFile, List<VCF> vcfList){
        System.out.println("LOG_INFO: writeVCFFile: execution started");
        log.log("LOG_INFO: writeVCFFile: execution started");
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter(pathToVCFFile, false));
            writer.write("@@format=VCF\n");
            writer.write("@@author=Kamran Abbaszade\n");
            writer.write("@@nickname=4LT4IR\n");
            writer.write("@@CHR\tPOS\tID\tREF\tALT\tQUALITY\tFORMAT\tMETADATA\n");


            for(VCF vcf: vcfList){
                writer.write(vcf.toString());
            }

            writer.close();
        }catch (IOException ioe){
            System.out.println(ioe.getMessage());
        }
        System.out.println("LOG_INFO: writeVCFFile: execution completed");
        log.log("LOG_INFO: writeVCFFile: execution completed");
    }

}



