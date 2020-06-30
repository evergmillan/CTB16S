/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Cut_region_804_1392;

import Alignment.Local_alignment;
import FastaIO.Fasta_Read;
import FastaIO.Fasta_Write;
import java.util.Iterator;
import java.util.LinkedHashMap;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 *
 * @author ever
 */
public class Primer_alignment {
    
    public String [] primerF804 = new String[6];
    public String [] primerR1392=new String[2];
    public int total=0,pequeñas=0,recortadas=0,noF=0,soloF=0;
    Local_alignment alineamiento=new Local_alignment();
    Fasta_Read ir=new Fasta_Read();
    Fasta_Write iw=new Fasta_Write();
    
    public void read_Map_seq(String file) throws CompoundNotFoundException
    {
        primerF804[0]="AGATTAGATACCCAAGTAGT";
        primerF804[1]="AGATTAGATACCCAGGTAGT";
        primerF804[2]="AGATTAGATACCCGAGTAGT";
        primerF804[3]="AGATTAGATACCCGGGTAGT";
        primerF804[4]="AGATTAGATACCCTAGTAGT";
        primerF804[5]="AGATTAGATACCCTGGTAGT";
        //Primers Reverse
        primerR1392[0]="ACGGGCGGTGTGTAC";
        primerR1392[1]="ACGGGCGGTGTGTGC";
        
        LinkedHashMap<String, String> a=ir.readFasta(file);
        LinkedHashMap<String, String> cut_seqs=new LinkedHashMap<>();
        Iterator<String> it = a.keySet().iterator();
            while(it.hasNext())
            {
                String id=it.next();
                String seq=a.get(id);
                String recortada=this.Sequence_cut(id,seq);
                if(recortada!=null)
                {
                    
                    String id_recortado=recortada.substring(0, recortada.indexOf("ñ"));
                    seq=recortada.substring(recortada.indexOf("ñ")+1);
                    if(seq.length()>550 && seq.length()<650)
                    {
                        if(cut_seqs.isEmpty())
                        {
                            System.out.println("Secuencia recortada no. "+recortadas +" id: "+id);
                            cut_seqs.put(seq, id_recortado);
                        }
                        if(!cut_seqs.containsKey(seq))
                        {
                            System.out.println("Secuencia recortada no. "+recortadas +" id: "+id);
                            cut_seqs.put(seq, id_recortado);
                        }
                    }                    
                }
                total++;
            }
            System.out.println("Mapa generado con "+cut_seqs.size()+" secuencias recortadas");
            Iterator<String> it_cut = cut_seqs.keySet().iterator();
            while(it_cut.hasNext())
            {
                String seq=it_cut.next();
                String id=cut_seqs.get(seq);
                
                iw.writeAst(id, seq, "smithella_802-1392");
            }
            
            System.out.println("Total: "+total);
            System.out.println("Recortadas: "+recortadas);
            System.out.println("Pequeñas: "+pequeñas);
            System.out.println("NOF: "+noF);
            System.out.println("SOLOF: "+soloF);
    }
    
    public String Sequence_cut(String id,String seq) throws CompoundNotFoundException
    {
        int pF,pR;
        pF=this.alignmentF(seq);
        if(pF>0)
        {
            String reverso=alineamiento.InversoComplementario(seq);
            pR=this.alignmentR(reverso);
            if(pR>pF)
            {
                String recortada=seq.substring(pF,pR+1);
                if(recortada.length()>550)
                {
                    recortadas++;
                    return id+" "+pF+" "+pR+" "+recortada.length()+"ñ"+recortada;
                }
                else
                {
                    pequeñas++;
                    return null;
                }
            }
            else
            {
                //System.out.println(reverso);
                soloF++;
            }
        }
        else
        {
            noF++;
        }
        return null;
    }
    
    public int alignmentF(String linea) throws CompoundNotFoundException
    {
        int rF=0;
        for(int j=0;j<primerF804.length;j++)
        {
            rF=alineamiento.extraerPosF(primerF804[j],linea);
            if(rF>0)
            {
                return rF;
            }
        }     
        return 0;
    }
    
    public int alignmentR(String linea) throws CompoundNotFoundException
    {
        int rR=0;
        for(int j=0;j<primerR1392.length;j++)
        {
            rR=alineamiento.extraerPosR(primerR1392[j],linea);
            if(rR>0)
            {
                return rR;
            }
        }     
        return 0;
    }
    
}
