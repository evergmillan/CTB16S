/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Alignment;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;

/**
 *
 * @author ever
 */
public class Local_alignment {
    
    public SequencePair<DNASequence, NucleotideCompound> alinearLocal(String querySeq,String targetSeq) throws CompoundNotFoundException
    {
        DNASequence target = new DNASequence(targetSeq,  
               AmbiguityDNACompoundSet.getDNACompoundSet());
       
        DNASequence query = new DNASequence(querySeq,  
               AmbiguityDNACompoundSet.getDNACompoundSet());

        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();  
         
        SimpleGapPenalty gapP = new SimpleGapPenalty();  
        gapP.setOpenPenalty((short)-10);  
        gapP.setExtensionPenalty((short)-4);  
         
        SequencePair<DNASequence, NucleotideCompound> psa =  
               Alignments.getPairwiseAlignment(query, target,  
                       Alignments.PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
        return psa;        
    }
    
    public String InversoComplementario(String linea) throws CompoundNotFoundException
    {
        linea=linea.replaceAll("[^ATGCN]","N");
        DNASequence seq=new DNASequence(linea);
        SequenceView reverso=seq.getReverseComplement();
        seq.setCompoundSet(reverso.getCompoundSet());
            
        return reverso.getSequenceAsString();
    }
    
    public int extraerPosF(String query, String target) throws CompoundNotFoundException
    {
        Local_alignment alinear=new Local_alignment();
        SequencePair<DNASequence, NucleotideCompound> psa=alinear.alinearLocal(query, target);
        if((((float)psa.getNumIdenticals()/(float)query.length())*100)>=85.0)
        {
            
            //Obtengo el primer alineado
            String primerA=psa.toString().substring(0,psa.toString().indexOf("\n"));
            //Obtengo la posicion de donde inicia el primer en la secuencia
            int posF=psa.getIndexInTargetForQueryAt(2);
            //Recorto el primer del alineamiento de la secuencia, menos 1 por que empieza en 0 y tamaño del primer alineado
            String primerC=target.substring(posF-1);
            if(primerA.length()>primerC.length())
            {
                primerC=primerC.substring(0,primerA.length()-1);
            }
            else
            {
                primerC=primerC.substring(0,primerA.length());
            }
            
            //Si termina en GT
            String ultimas2=primerC.substring(primerC.length()-2);
            if(ultimas2.compareTo("GT")==0)
            {
                return posF-1;
            }
        }
        return 0;
    }
    
    public int extraerPosR(String query, String target) throws CompoundNotFoundException
    {
        String reversoCP=query;
        String primerC="";
        SequencePair<DNASequence, NucleotideCompound> psa=this.alinearLocal(reversoCP, target);
        if((((psa.getNumIdenticals()/(float)query.length())*100)>=93.0))
        {
            //Obtengo el primer alineado
            String primerA=psa.toString().substring(0,psa.toString().indexOf("\n"));
            
            //Obtengo la posicion de donde inicia el primer en la secuencia
            int posR=psa.getIndexInTargetForQueryAt(2);
            
            if(posR==1)
            {
                //Recorto el primer del alineamiento de la secuencia, menos 1 por que empieza en 0 y tamaño del primer alineado
                primerC=target.substring((posR-1),(posR-1)+primerA.length());
            }
            else
            {
                //Recorto el primer del alineamiento de la secuencia, menos 1 por que empieza en 0 y tamaño del primer alineado
                primerC=target.substring((posR-2),(posR-2)+primerA.length());
            }
            
            //primerC=primerC.substring(0,primerA.length()-1);
            
            //Si termina en AC || GC
            String ultimas2=primerC.substring(primerC.length()-2);
            
            if(ultimas2.compareTo("AC")==0 || ultimas2.compareTo("GC")==0)
            {
                //Retorno la posicion donde termina el primer reverso
                return ((target.length()-1-posR)-(primerA.length()-1));
            }
        }
        return 0;
    }
    
}
