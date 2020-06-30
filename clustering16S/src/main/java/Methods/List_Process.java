/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Methods;

import FastaIO.Fasta_Read;
import FastaIO.Fasta_Write;
import Objects_16S.Taxo_Seq;
import java.util.ArrayList;
import java.util.Comparator;
    import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
 
/**
 *
 * @author ever
 */
public class List_Process {
    
    @SuppressWarnings("UnusedAssignment")
    public Double[] convertir_seq_Double (String seq)
    {
        Double[] seq_dtw= new  Double[seq.length()];
        for (int i = 0; i < seq.length(); i++) 
        {
            if(seq.charAt(i)=='A')
            {
               seq_dtw[i]=1.0;
            }
            if(seq.charAt(i)=='T')
            {
               seq_dtw[i]=2.0;
            }
            if(seq.charAt(i)=='G')
            {
               seq_dtw[i]=3.0;
            }
            if(seq.charAt(i)=='C')
            {
               seq_dtw[i]=4.0;
            }
            if(seq.charAt(i)=='N')
            {
               seq_dtw[i]=5.0;
            }            
        }
        return seq_dtw;
    }
    public Double[] convertir_protein_Double (String seq)
    {
        Double[] seq_dtw= new  Double[seq.length()];
        for (int i = 0; i < seq.length(); i++) 
        {
            if(seq.charAt(i)=='A')
            {
               seq_dtw[i]=1.0;
            }
            if(seq.charAt(i)=='R')
            {
               seq_dtw[i]=2.0;
            }
            if(seq.charAt(i)=='N')
            {
               seq_dtw[i]=3.0;
            }
            if(seq.charAt(i)=='D')
            {
               seq_dtw[i]=4.0;
            }
            if(seq.charAt(i)=='C')
            {
               seq_dtw[i]=5.0;
            }
            if(seq.charAt(i)=='Q')
            {
               seq_dtw[i]=6.0;
            }
            if(seq.charAt(i)=='E')
            {
               seq_dtw[i]=7.0;
            }
            if(seq.charAt(i)=='G')
            {
               seq_dtw[i]=8.0;
            }
            if(seq.charAt(i)=='H')
            {
               seq_dtw[i]=9.0;
            }
            if(seq.charAt(i)=='I')
            {
               seq_dtw[i]=10.0;
            }
            if(seq.charAt(i)=='L')
            {
               seq_dtw[i]=11.0;
            }
            if(seq.charAt(i)=='K')
            {
               seq_dtw[i]=12.0;
            }
            if(seq.charAt(i)=='M')
            {
               seq_dtw[i]=13.0;
            }
            if(seq.charAt(i)=='F')
            {
               seq_dtw[i]=14.0;
            }
            if(seq.charAt(i)=='P')
            {
               seq_dtw[i]=15.0;
            }
            if(seq.charAt(i)=='S')
            {
               seq_dtw[i]=16.0;
            }
            if(seq.charAt(i)=='T')
            {
               seq_dtw[i]=17.0;
            }
            if(seq.charAt(i)=='W')
            {
               seq_dtw[i]=18.0;
            }
            if(seq.charAt(i)=='Y')
            {
               seq_dtw[i]=19.0;
            }
            if(seq.charAt(i)=='V')
            {
               seq_dtw[i]=20.0;
            }
            
        }
        
        //String seq2=seq.replaceAll("[^ARNDCQEGHILKMFPSTWYV]", "ñ");
        return seq_dtw;
    }
    public Double distanciaHammit(double a,double b)
    {
        Double distancia=0.0;
        if(a==b)
        {
            return distancia;
        }
        else
        {
          return 1.0;  
        }
    }
    
    @SuppressWarnings({"static-access", "SizeReplaceableByIsEmpty"})
    public void verify_seq_to_annotation() throws CompoundNotFoundException
    {
        Fasta_Read read=new Fasta_Read();
        List<List<Taxo_Seq>> a=read.readFastaListBDcode("BDcode_v3_taxa.fasta");
        List<Taxo_Seq> BD_code_v3_final=new ArrayList();
        Fasta_Write write=new Fasta_Write();
        /*a.sort(Comparator.comparing(Taxo_Seq::getSequence));
        for (int i = 0; i < a.size(); i++) {
            System.out.println(a.get(i).getGenre());
        }*/
        int ban=0;
        
        @SuppressWarnings("UnusedAssignment")
        //Iterador para recorrer la lista de listas (por género)
        Iterator <List<Taxo_Seq>>  itr = a.iterator(); 
        while (itr.hasNext())
        {
            //Se crea un objeto para guardar el elemento actual
            @SuppressWarnings("UnusedAssignment")
            Taxo_Seq now=new Taxo_Seq();
            @SuppressWarnings("UnusedAssignment")
            Taxo_Seq aux=new Taxo_Seq();
            List<Taxo_Seq> group=itr.next();
            itr.remove();
            int pos=0;
            //Se calculan las distncias todos contra todos de cada género, se guardan aquellas que estan entre el 97 y 100
            while (pos+1<group.size())
            {
                for (int i = pos+1; i < group.size(); i++) 
                {
                    now=group.get(pos);
                    aux=group.get(i);
                    if(this.distance_calculate(now, aux))
                    {
                        ban=this.exists_seqs(BD_code_v3_final, now, aux);
                        if(ban==2)
                        {
                            BD_code_v3_final.add(now);
                            BD_code_v3_final.add(aux);
                            System.out.println(now.getTaxa());
                            System.out.println(aux.getTaxa());
                        }
                        
                        if(ban==10)
                        {
                            BD_code_v3_final.add(now);
                            System.out.println(now.getTaxa());
                        }
                        
                        if(ban==20)
                        {
                            BD_code_v3_final.add(aux);
                            System.out.println(aux.getTaxa());
                        }
                    }
                }
                pos++;
            }
            
            ban=0;
        }
        String id="",sequence="";
        for (int i = 0; i < BD_code_v3_final.size(); i++) 
        {
            id=BD_code_v3_final.get(i).getId()+" "+BD_code_v3_final.get(i).getTaxa();
            sequence=BD_code_v3_final.get(i).getSequence();
            write.writeAst(id, sequence, "BD_code_v3_final");
        }
        System.out.println(BD_code_v3_final.size());
    }
    
    @SuppressWarnings("static-access")
    public boolean distance_calculate(Taxo_Seq now,Taxo_Seq after)
    {
        DTW dtw=new DTW();
        double distance=0.0;
        Double dwr[][] = new Double[this.convertir_seq_Double(now.getSequence()).length][this.convertir_seq_Double(after.getSequence()).length];
        Stack<Double> wr = new Stack<>();
        distance=dtw.dtw(this.convertir_seq_Double(now.getSequence()),this.convertir_seq_Double(after.getSequence()),dwr,wr, " ");
        if(distance>97 && distance <=100)
        {
            return true;
        }
        return false;
    }
    
    public int exists_seqs(List<Taxo_Seq> data, Taxo_Seq now, Taxo_Seq after)
    {
        int now_int=10,after_int=20,ban_now=0,ban_after=0;
        if(data.isEmpty())
        {
            //Si se retorna 2 es porque no existen los dos elementos o esta vacía la lista
            return 2;
        }
        else
        {
            for (int i = 0; i < data.size(); i++) 
            {
                //Se evalua si existe el elemento actual
                if(data.get(i).getId().equals(now.getId()))
                {
                    //Si existe el elemento actual, ban_now cambia a 1
                    ban_now=1;
                }
                //Se evalua si existe el elemento siguiente
                if(data.get(i).getId().equals(after.getId()))
                {
                    ////Si existe el elemento siguiente, ban_after cambia a 1
                    ban_after=1;
                    //Si existe el elemento actual, ademas del elemento anterior
                    if(ban_now==1)
                    {
                        //Se para el ciclo, ya que existen los dos elementos
                        return 0;
                    }                    
                }
            }
        }
        
        //Se evaluan los elementos que existen en la lista
        if(ban_now==0 && ban_after==0)
        {
            return 2;
        }
        if(ban_now==1 && ban_after==0)
        {
            return after_int;
        }
        
        if(ban_now==0 && ban_after==1)
        {
            return now_int;
        }
        
        return 0;
    }
    
    
    //Funcion para agrupar datos contra los centroides formados por la BD RefSeq
    //recibe la lista con los centroides y semillas, y como segundo parámetro la lista de secuencias a agrupar
    @SuppressWarnings({"static-access", "UnusedAssignment"})
    public void clusteringData(List<List<Taxo_Seq>> centroids, List<Taxo_Seq> data)
    {
        //Lista de listas con las distancias para cada secuencia contra cada centroide
        List<List<Double>> averages_calculates=new ArrayList();
        //Lista con los promedios para cada grupo con las secuencias asignadas en cada iteración
        List<Double> averages_groups=new ArrayList();
        //Objeto con los metodos para medir distancia entre dos secuencias
        DTW dtw=new DTW();
        //Lista para guardar las distancias para cada secuencia que se agrupará
        @SuppressWarnings("UnusedAssignment")
        double distance = 0.0;
        //Variable para guardar la posición del genero con mayor similitud con respecto a la lista de centroides
        int pos=0;
        //Variable para guardar la distancia entre dos secuencias
        double promedio=0.0,tam;
        //Pila para guardar los valores de los metodos de DTW
        Stack<Double> wr = new Stack<>();
        
        //---Inicia Kmeans---//
        
        //---Primera iteración, asignación de grupo por primera vez---//
        
        //Ciclo para recorrer las secuencias a agrupar
        for (int i = 0; i < data.size(); i++) 
        {
            //Se crea la lista para las distancias promedio entre la secuencia de entrada y las semillas de los centroides
            @SuppressWarnings("MismatchedQueryAndUpdateOfCollection")
            List<Double> averages=new ArrayList();
            
            //Ciclo para recorrer la lista de listas de los centroides
            for (int j = 0; j < centroids.size(); j++) 
            {
                //Se obtiene la cantidad de semillas para cada centroid
                tam=centroids.get(j).size();
                
                //Ciclo para recorrer las semillas por lista
                for (int k = 0; k < centroids.get(j).size(); k++) 
                {
                    //----Proteina----
                    /*
                    //Se crea la matriz con los diferentes valores para cada base nitrogenada (A-T-G-C)
                    Double dwr[][] = new Double[this.convertir_protein_Double(data.get(i).getProtein_seq().substring(0,100)).length][this.convertir_protein_Double(centroids.get(j).get(k).getProtein_seq().substring(0,100)).length];
                    //Se obtiene la distancia entre secuencia de entrada y cada semilla
                    distance=dtw.dtw(this.convertir_protein_Double(data.get(i).getProtein_seq().substring(0,100)),this.convertir_protein_Double(centroids.get(j).get(k).getProtein_seq().substring(0,100)),dwr,wr);
                    //System.out.println(data.get(i).getId()+"\t"+centroids.get(j).get(k).getId()+"\t"+distance);
                    //Se suman las distancias y se guardan en la variable promedio
                    promedio=promedio+distance;*/
                    
                    //----Nucleotido----
                    //Se crea la matriz con los diferentes valores para cada base nitrogenada (A-T-G-C)
                    Double dwr[][] = new Double[this.convertir_seq_Double(data.get(i).getSequence()).length][this.convertir_seq_Double(centroids.get(j).get(k).getSequence()).length];                    
                    //Se obtiene la distancia entre secuencia de entrada y cada semilla
                    distance=dtw.dtw(this.convertir_seq_Double(data.get(i).getSequence()),this.convertir_seq_Double(centroids.get(j).get(k).getSequence()),dwr,wr,centroids.get(j).get(k).getGenre());
                    //System.out.println(data.get(i).getId()+"\t"+centroids.get(j).get(k).getId()+"\t"+distance);
                    //Se suman las distancias y se guardan en la variable promedio
                    dwr=null;
                    dtw=null;
                    promedio=promedio+distance;
                    if(distance==0.0)
                    {
                        k=centroids.get(j).size();
                    }
                }
                
                //Se guarda el promedio de distancias entre la secuencia data_i y centroide j
                promedio=promedio/tam;
                if(promedio>98)
                {
                    averages.add(promedio);
                    //Se reinicializa el promedio
                    promedio=0;
                    break;
                }
                else
                {
                    /*if(promedio>80)
                    {
                        averages.add(promedio);
                    }*/
                    averages.add(promedio);
                    //Se reinicializa el promedio
                    promedio=0;
                }
                
            }
                      
            //Se extrae el genero del centroide que tiene mas similitud con la secuencia de entrada y se asigna a la secuencia de entrada
            pos=this.extract_max_identity(averages);
            data.get(i).setGrupo_genre(centroids.get(pos).get(0).getGenre());
            data.get(i).setPos_genre(pos);
            data.get(i).setDistance(averages.get(pos));
            if(data.get(i).getId().contains(data.get(i).getGrupo_genre()))
            {
                System.out.println(data.get(i).getDistance()+"\t"+data.get(i).getGrupo_genre()+"\t"+data.get(i).getId());
            }
            else
            {
                System.out.println("");
                System.out.println("------------------");
                System.out.println(data.get(i).getDistance()+"\t"+data.get(i).getGrupo_genre()+"\t"+data.get(i).getId());
                System.out.println("------------------");
                System.out.println("");
            }
            
            
            //Se agrega la lista de distancias calculadas a la lista de listas
            averages_calculates.add(averages);
            averages.clear();
            averages=null;
            Runtime garbage = Runtime.getRuntime();
            garbage.gc();
        }
        
        //---Termina primera iteración---//
        
        //Variable para saber si hay algun cambio
        int cambio=0, iteraciones=0;
        
        //Se ordena por genero las secuencias ya asignadas a un grupo
        data.sort(Comparator.comparing(Taxo_Seq::getPos_genre));
        
        //---Inicia el recalculo de distancias para cada grupo (recalculo de centroides)---//
        
        //Hacer mientras haya cambios de grupos
       /* do
        {
            
            //averages_groups
            iteraciones++;
        }while(cambio!=0);*/
        
        
        
        
    }
  
    
    public List<Double> average_calculate(List<Taxo_Seq> data, int num_centroids)
    {
        //Lista para guardar los promedios por grupo
        List<Double> averages_groups=new ArrayList();
        
        double average=0.0, aux=0.0,suma=0.0;
        int limite=0,cont=0;
        //Para cada centroide calcular el promedio de las secuencias asignadas
        for (int i = 0; i < num_centroids; i++) 
        {
            //Ciclo para recorrer la lista de objetos del tipo Taxo_Seq
            for (int j = limite; j < data.size(); j++) 
            {
                //Si son del mismo grupo, se empieza en 1 ya que se recibe la lista data ordenada por numero de grupo
                if(data.get(j).getPos_genre()==i)
                {
                    //Se suma la distance para las secuencias del mismo grupo
                    suma=suma+data.get(j).getDistance();
                    
                    //Se incrementa el limite del grupo
                    limite++;
                    //Se incrementa el numero de miembros en el grupo
                    cont++;
                }
                else
                {
                    //La secuencia pertenece a otro grupo, se sale del ciclo
                    j=data.size();
                }
            }
            
            //Se calcula el promedio por grupo
            average=(suma/(double)cont);
            //Se agrega el promedio a la lista de promedios
            averages_groups.add(average);
            //Se reinicializan las variables
            cont=0;
            suma=0;
            
        }
        return averages_groups;
    }
    
    public int extract_max_identity(List<Double> averages)
    {
        double max_num=0.0;
        int pos=0;
        for(int i=0; i<averages.size(); i++)
        {            
            if(averages.get(i)>max_num)
            {
                max_num = averages.get(i);
                pos=i;
            }
        }        
        return pos;
    }
    
    public String translate_frame_1(String seq) throws CompoundNotFoundException
    {
        DNASequence dna = new DNASequence(seq);
        TranscriptionEngine te = TranscriptionEngine.getDefault();
        Frame[] frames = Frame.getForwardFrames();
        Map<Frame, Sequence<AminoAcidCompound>> results = te.multipleFrameTranslation(dna, frames);
        Iterator it = results.keySet().iterator();
        String protein_frame_1="";
        int cont=1;
        while(it.hasNext()){
            Frame key = (Frame) it.next();
            protein_frame_1=results.get(key).getSequenceAsString().replaceAll("\\*", "");
            //return protein_frame_1.replaceAll("X", "");
            if(cont==1)
            {
                return protein_frame_1.replaceAll("X", "");
            }
            cont++;            
        }
        return null;
    }
}
