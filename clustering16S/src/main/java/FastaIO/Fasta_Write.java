/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaIO;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author cisei31
 */
public class Fasta_Write {
    
    public void writeAst(String id, String secuencia,String nombre) {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(nombre+".fasta",true));)
        { 
            //Escribimos en el fichero
            bw.write(">"+id);
            bw.newLine();
            bw.write(secuencia);
            bw.newLine();
            //Guardamos los cambios del fichero
            bw.flush();
                      
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
    }
    
    public void writeSeeds(String linea,String nombre) {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(nombre+".txt",true));)
        { 
            //Escribimos en el fichero
            bw.write(linea);
            bw.newLine();
            //Guardamos los cambios del fichero
            bw.flush();
                      
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
    }
    
}
