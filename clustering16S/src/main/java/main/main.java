/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import Cut_region_804_1392.Primer_alignment;
import FastaIO.Fasta_Read;
import Methods.List_Process;
import Objects_16S.Taxo_Seq;
import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 *
 * @author cisei31
 */
public class main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws CompoundNotFoundException {
        
        //Primer_alignment align_primer=new Primer_alignment();
        //align_primer.read_Map_seq("sequence.fasta");
        List_Process process=new List_Process();
        
        //process.extract_family_list();
        //process.verify_seq_to_annotation();
        Fasta_Read read=new Fasta_Read();
        List<List<Taxo_Seq>> b=read.readFastacentroidSeeds("centroids_97-100.fasta");
        List<Taxo_Seq> a=read.readFastaList("prueba.fasta");
        //List<Taxo_Seq> a=read.readFastaList("BD_804-1392_refSeq.fasta");
        process.clusteringData(b, a);
       // read.readFastaList("BDcode_v2.fasta");
    }
    
}
