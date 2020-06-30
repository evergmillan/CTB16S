/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FileIO;

import FastaIO.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
/**
 *
 * @author morelos
 */
public class File_Read {
    
            
    public final String url="/home/ever/Documentos/Clasificados/BD_Control/Analisis_Consenso/Pruebas/";
    
    public LinkedHashMap<String, String> readFile(String file) throws CompoundNotFoundException 
    {
        LinkedHashMap<String, String> a = new LinkedHashMap<>();
        try(BufferedReader br=new BufferedReader(new FileReader( file ));)
        {
            String id = null;
            String linea;
            String seq ="",anterior="";
            while(br.ready())
            {
                linea=br.readLine();
                if(linea.startsWith(">"))
                {
                    id=linea;
                    if(seq.compareTo("")!=0)
                    {
                        a.put(anterior, seq);
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
                a.put(id, seq);
            }            
            System.out.println(a.size());
            return a;
            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return null;
    }
}