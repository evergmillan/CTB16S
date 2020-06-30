/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaIO;

import Methods.List_Process;
import Objects_16S.Taxo_Seq;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
/**
 *
 * @author morelos
 */
public class Fasta_Read {
    private String ID;
    private String seq;

    public Fasta_Read(String ID, String seq) {
        this.ID = ID;
        this.seq = seq;
    }

    public Fasta_Read() {
    }
    
    public String getID() {
        return ID;
    }

    public void setID(String ID) {
        this.ID = ID;
    }

    public String getSeq() {
        return seq;
    }

    public void setSeq(String seq) {
        this.seq = seq;
    }
    
    
            
    public final String url="/home/ever/Documentos/Clasificados/BD_Control/Analisis_Consenso/Pruebas/";
    
    public LinkedHashMap<String, String> readFasta(String file) throws CompoundNotFoundException 
    {
        LinkedHashMap<String, String> a = new LinkedHashMap<>();
        try(BufferedReader br=new BufferedReader(new FileReader( file ));)
        {
            String id = null;
            String linea,taxo;
            @SuppressWarnings("LocalVariableHidesMemberVariable")
                    String seq = "";
            String anterior = "";
            while(br.ready())
            {
                linea=br.readLine();
                if(linea.startsWith(">"))
                {
                    id=linea;
                    if(seq.compareTo("")!=0)
                    {
                        taxo=anterior.substring(1,anterior.indexOf(" ")+1);
                        //String [] id_full=taxo.split("|");
                        //if(id_full.length==7)
                        //{
                            a.put(taxo, seq.replaceAll("U", "T"));                        
                        //}
                        seq="";
                        
                    }            
                }
                else
                {
                    seq+=linea.substring(0, linea.length());
                    anterior=id;                    
                }
            }
            if(br.ready()==false)
            {
                taxo=anterior.substring(1,anterior.indexOf(" ")+1);
                //String [] id_full=taxo.split(";");
                //if(id_full.length==7)
                //{
                    a.put(taxo, seq.replaceAll("U", "T"));                                
                //}
            }
            return a;
            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return null;
    }
    
    @SuppressWarnings("UnusedAssignment")
    public List<Taxo_Seq> readFastaList(String file) throws CompoundNotFoundException 
    {
        List_Process process=new List_Process();
        @SuppressWarnings("UnusedAssignment")
        List<Taxo_Seq> a = new ArrayList<>();
        String family ="",genre="",id_taxa="",orden="";
        
        String description="";
        
        
        try(BufferedReader br=new BufferedReader(new FileReader( file ));)
        {
            @SuppressWarnings("UnusedAssignment")
            String [] id;
            String linea;
            @SuppressWarnings("LocalVariableHidesMemberVariable")
            String seq = "";
            while(br.ready())
            {
                linea=br.readLine();
                if(linea.startsWith(">"))
                {
                    if(seq.length()>0)
                    {
                        Taxo_Seq element=new Taxo_Seq();
                        element.setId(id_taxa);
                        element.setSequence(seq.replaceAll("[^ATGCN]", "N"));
                        element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                        a.add(element);
                        
                        seq="";
                    }
                    linea=linea.substring(1);                    
                    id_taxa=linea;
                }
                else
                {
                    seq+=linea.substring(0, linea.length());                    
                }
            }
            //Se agrega la ultima secuencia a la lista de listas
            if(br.ready()==false)
            {
                Taxo_Seq element=new Taxo_Seq();
                element.setId(id_taxa);
                element.setSequence(seq.replaceAll("[^ATGCN]", "N"));
                element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                a.add(element);
                seq="";
            }            
            System.out.println("---Load data finish---");
            return a;
            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return null;
    }
    
    public List<List<Taxo_Seq>> readFastaListBDcode(String file) throws CompoundNotFoundException 
    {
        List_Process process=new List_Process();
        @SuppressWarnings("UnusedAssignment")
        List<Taxo_Seq> a = new ArrayList<>();
        String genre="";       
        
        try(BufferedReader br=new BufferedReader(new FileReader( file ));)
        {
            @SuppressWarnings("UnusedAssignment")
            String id = null,taxon = null;
            String [] taxa;
            String linea;
            @SuppressWarnings("LocalVariableHidesMemberVariable")
            String seq = "";
            String regex = "^[A-Z][A-Z]";
            Pattern patron = Pattern.compile(regex);
            while(br.ready())
            {
                linea=br.readLine();
                if(linea.startsWith(">"))
                {
                    if(seq.length()>0)
                    {
                        
                        Taxo_Seq element=new Taxo_Seq();
                        element.setGenre(genre);
                        element.setId(id);
                        element.setSequence(seq.replaceAll("[^ATGCN]", "N"));
                        element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                        element.setTaxa(taxon);
                        Matcher emparejador = patron.matcher(element.getGenre());
                        boolean esCoincidente = emparejador.find();
                        if (!esCoincidente) 
                        {
                            a.add(element);
                        }
                        
                        
                        seq="";
                    }
                    linea=linea.substring(1);
                    taxon=linea.substring(linea.indexOf(" ")+1);
                    id=linea.substring(0,linea.indexOf(" "));
                    linea=linea.substring(linea.indexOf(" ")+1);
                    taxa=linea.split(",");
                    //orden=id[3];
                    //family=id[2];
                    genre=taxa[1];
                }
                else
                {
                    seq+=linea.substring(0, linea.length());                    
                }
            }
            //Se agrega la ultima secuencia a la lista de listas
            if(br.ready()==false)
            {
                Taxo_Seq element=new Taxo_Seq();
                //element.setOrden(orden);
                //element.setFamily(family);
                element.setGenre(genre);
                element.setId(id);
                element.setSequence(seq);
                element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                element.setTaxa(taxon);
                Matcher emparejador = patron.matcher(element.getGenre());
                boolean esCoincidente = emparejador.find();
                if (!esCoincidente) 
                {
                    a.add(element);
                }
                seq="";
            }
            //Se crea la lista de listas por genero
            List<List<Taxo_Seq>> genres=new ArrayList<>();
            int ban=0;
            Iterator<Taxo_Seq> it = a.iterator();
            while (it.hasNext()) 
            {
                if(genres.isEmpty())
                {
                     @SuppressWarnings("UnusedAssignment")
                    Taxo_Seq element=new Taxo_Seq();
                    element=it.next();                   
                    List<Taxo_Seq> group=new ArrayList<>();
                    group.add(element);
                    genres.add(group);
                    //it.remove();
                }
                else
                {
                     @SuppressWarnings("UnusedAssignment")
                    Taxo_Seq element=new Taxo_Seq();
                    element=it.next();
                    for (int i = 0; i < genres.size(); i++) 
                    {
                        if(genres.get(i).get(0).getGenre().equals(element.getGenre()))
                        {
                            genres.get(i).add(element);
                            //it.remove();
                            ban=1;
                            break;
                        }
                    }
                    
                    if(ban==0)
                    {
                        List<Taxo_Seq> group=new ArrayList<>();
                        group.add(element);
                        genres.add(group);
                        //it.remove();
                    }
                    ban=0;
                }
            }
            Iterator <List<Taxo_Seq>>  itr = genres.iterator(); 
            while (itr.hasNext()) 
            { 
                if(itr.next().size()<2)
                {
                    itr.remove();
                }
            }
            System.out.println("---Load data finish---");
            return genres;
            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return null;
    }
    
    @SuppressWarnings("ImplicitArrayToString")
    public List<List<Taxo_Seq>> readFastacentroidSeeds(String file) throws CompoundNotFoundException 
    {        
        @SuppressWarnings("UnusedAssignment")
        List<List<Taxo_Seq>> centroids = new ArrayList<>();
        List<Taxo_Seq> a = new ArrayList<>();
        String genre="",id_taxa="";
        int ban=0;
        List_Process process=new List_Process();
        
        try(BufferedReader br=new BufferedReader(new FileReader( file ));)
        {
            @SuppressWarnings("UnusedAssignment")
            String [] id_array;
            String linea,id,taxa="";
            @SuppressWarnings("LocalVariableHidesMemberVariable")
            String seq = "";
            while(br.ready())
            {
                linea=br.readLine();
                if(linea.startsWith(">"))
                {
                    if(seq.length()>0)
                    {
                        Taxo_Seq element=new Taxo_Seq();

                        element.setGenre(genre);
                        element.setId(id_taxa);
                        element.setTaxa(taxa);
                        element.setSequence(seq.replaceAll("[^ATGCN]", "N"));
                        element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                        if(a.isEmpty())
                        {
                            a.add(element);
                            if(centroids.isEmpty())
                            {
                                centroids.add(a);
                            }
                        }
                        else
                        {
                            for (int i = 0; i < centroids.size(); i++) 
                            {
                                if(element.getGenre().equals(centroids.get(i).get(0).getGenre()))
                                {                             
                                    centroids.get(i).add(element);
                                    ban=1;
                                    break;
                                }
                            }
                            
                            if(ban==0)
                            {
                                a = new ArrayList<>();
                                a.add(element);                                
                                centroids.add(a);
                            }
                        }
                        ban=0;
                        
                        seq="";
                    }
                    linea=linea.substring(1);
                    taxa=linea.substring(linea.indexOf(" ")+1);
                    id_taxa=linea.substring(0, linea.indexOf(" "));
                    linea=linea.substring(linea.indexOf(" ")+1);
                    id_array=linea.split(",");
                    genre=id_array[1];
                    //id_taxa=id_array[0].substring(0, id_array[0].indexOf(" "));
                    
                }
                else
                {
                    seq+=linea.substring(0, linea.length());                    
                }
            }
            ban=0;
            //Se agrega la ultima secuencia a la lista de listas
            if(br.ready()==false)
            {
                Taxo_Seq element=new Taxo_Seq();
                element.setGenre(genre);
                element.setId(id_taxa);
                element.setTaxa(taxa);
                element.setSequence(seq);
                element.setProtein_seq(process.translate_frame_1(element.getSequence()));
                for (int i = 0; i < centroids.size(); i++) 
                {
                    if(element.getGenre().equals(centroids.get(i).get(0).getGenre()))
                    {              
                        centroids.get(i).add(element);                        
                        ban=1;
                        break;
                    }
                }
                if(ban==0)
                {
                    a = new ArrayList<>();
                    a.add(element);                                
                    centroids.add(a);
                }
                seq="";
            }
            System.out.println("---Load centroids finish---");
            return centroids;
            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return null;
    }
}