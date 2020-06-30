/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cut_region_804_1392;
import java.util.LinkedHashMap;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;

/**
 *
 * @author morelos
 */
public class Cut_region {
    
    public int alignment (String query, String primer, int tipo) throws CompoundNotFoundException
    {
        if(tipo==1)
        {
            GapPenalty penalty = new SimpleGapPenalty( -10, -4);
            PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(
            new DNASequence(query, AmbiguityDNACompoundSet.getDNACompoundSet()),
            new DNASequence(primer, AmbiguityDNACompoundSet.getDNACompoundSet()),
            PairwiseSequenceAlignerType.LOCAL,
            penalty, SubstitutionMatrixHelper.getNuc4_4());
            SequencePair<DNASequence, NucleotideCompound> alignment = aligner.getPair();
            int identical = alignment.getNumIdenticals();
            int pos1=alignment.toString().indexOf("\n");
            String primera=alignment.toString().substring(pos1-2,pos1);
            String ultima=alignment.toString().substring(alignment.toString().length()-3, alignment.toString().length()-1);
            
            if(((ultima.compareTo("GT")==0 && primera.compareTo("GT")==0) /*|| (primera.compareTo("GT")==0 && ultima.compareTo("GT")==0)*/) && identical / (float) primer.length()>=0.85)
            {               

                int posF=((alignment.getIndexInQueryForTargetAt(2)-2));
                if(posF<0)
                {
                    return -1;
                }

                return posF;                
            }
        }
        else
        {
            DNASequence seq=new DNASequence(query);
            SequenceView reverso=seq.getReverseComplement();
            seq.setCompoundSet(reverso.getCompoundSet());
            LinkedHashMap <String,AmbiguityDNACompoundSet> seqs=new LinkedHashMap<>();
            seqs.put(query, AmbiguityDNACompoundSet.getDNACompoundSet());
            GapPenalty penalty = new SimpleGapPenalty( -20, -5);
            PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(
            new DNASequence(reverso.getSequenceAsString(), AmbiguityDNACompoundSet.getDNACompoundSet()),
            new DNASequence(primer, AmbiguityDNACompoundSet.getDNACompoundSet()),
            PairwiseSequenceAlignerType.LOCAL,
            penalty, SubstitutionMatrixHelper.getNuc4_4());
            SequencePair<DNASequence, NucleotideCompound> alignment = aligner.getPair();
            int identical = alignment.getNumIdenticals();
            int pos1=alignment.toString().indexOf("\n");
            String primera = null;
            if(pos1-2>0)
            {
              primera=alignment.toString().substring(pos1-2,pos1);   
            }
            else
            {
                return -1;
            }
            
            String ultima=alignment.toString().substring(alignment.toString().length()-3, alignment.toString().length()-1);
            
            if((ultima.compareTo("AC")==0 && primera.compareTo("AC")==0) || (primera.compareTo("GC")==0 && ultima.compareTo("GC")==0))
            {
                if(identical / (float) primer.length()>=0.93)
                {
                    int posR=((query.length()-alignment.getIndexInQueryForTargetAt(2))+2);
                  
                    if(posR>query.length())
                    {
                        posR=query.length();
                    }
                    else
                    {
                        if(posR<0 || posR==0)
                            return -1;
                    }
                    return posR;
                    
                }               
            }            
        }
        return 0;
    }   
}
