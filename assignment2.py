#! /usr/bin/env python3
# vim: set fileencoding=UTF-8 :

import vcf

__author__ = 'mgollobich'


class Assignment2:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION, "\n")


        self.filename = 'AmpliseqExome.20141120.NA24385.vcf'
        self.file = vcf.Reader(open(self.filename, 'r'))
        #rec = next(self.file)
        
        #Used Reference: https://github.com/jamescasbon/PyVCF/blob/master/README.rst
        


    def get_average_quality_of_son(self):
        '''
        Get the average PHRED quality of all variants
        :return: 
        '''    
        print("Calculating avg record.QUAL")
        self.file = vcf.Reader(open(self.filename, 'r'))
        quality = 0
        n = 0
        for record in self.file:
            quality += record.QUAL
            n += 1
        avg = quality/n
        return avg, n

        
        
    def get_total_number_of_variants_of_son(self):
        '''
        Get the total number of variants
        :return: total number of variants
        '''
        number = 0
        for record in self.file:
            number += 1
        return number
    
    
    def get_variant_caller_of_vcf(self):
        '''
        Return the variant caller name
        :return: 
        '''
        vcaller = self.file.metadata['source'][0]
        return vcaller
        
        
    def get_human_reference_version(self):
        '''
        Return the genome reference version
        :return: 
        '''
        self.file = vcf.Reader(open(self.filename, 'r'))
        ## Jede Zeile des Headers (_header_lines) wird durchsucht
        for line in self.file._header_lines:
            ## wenn die Zeile mit ##reference beginnt, enthält sie die Info über den Caller
            ## aber da das öfter vorkommt, spezifiziere ich es zusätzlich
            ## das Format: ##reference=file:///results/referenceLibrary/tmap-f3/hg19/hg19.fasta
            if line.startswith('##reference=file://'):
                ## zuerst wird die zweite Hälfte entfernt indem man die Zeile an "." teilt
                version = line.split(".")
                ## damit man nur hg19 bekommt, splite ich den String nochmal
                version = version[0].split("/")
                ## beim Split wird eine Liste zurückgegeben deshalb verwende ich den 0. (1.) Eintrag.
        return version[7]
        
        
    def get_number_of_indels(self):
        '''
        Return the number of identified INDELs 
        :return: 
        '''
        reader = vcf.Reader(open(self.filename, 'r'))
        indel = 0
        for record in reader:
            if record.is_indel:
                indel += 1
        return indel

    def get_number_of_snvs(self):
        '''
        Return the number of SNVs
        :return: 
        '''
        reader = vcf.Reader(open(self.filename, 'r'))
        snp = 0
        for record in reader:
            if record.is_snp:
                snp += 1
        return snp
    def get_number_of_heterozygous_variants(self):
        '''
        Return the number of heterozygous variants
        :return: 
        '''
        reader = vcf.Reader(open(self.filename, 'r'))
        n_het = 0
        for record in reader:
            if record.num_het:
                n_het += 1
        return n_het
        
    def print_summary(self):
        print(self.get_average_quality_of_son()," - Average quality of SON, Number of Variants ")
        print("The variant caller of vcf is: ",self.get_variant_caller_of_vcf(),"\n")
        print('{:<8}'.format(self.get_number_of_indels()), " - Number of Indels")
        print('{:<8}'.format(self.get_number_of_snvs()), " - Number of SNPs")
        print('{:<8}'.format(self.get_number_of_heterozygous_variants())," - Number of Heterozygous Variants")   
        print('{:<8}'.format(self.get_human_reference_version())," - Version of reference genome")
        
        
if __name__ == '__main__':
    print("Assignment 2: ",__author__)
    assignment1 = Assignment2()
    assignment1.print_summary()

#Used Reference: https://github.com/jamescasbon/PyVCF/blob/master/README.rst

