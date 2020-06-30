/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Methods;

import java.util.Stack;

/**
 *
 * @author ever
 */
public class DTW {
    @SuppressWarnings("SuspiciousIndentAfterControlStatement")
    public static Double dtw(Double a[],Double b[],Double dw[][], Stack<Double> w, String genero)
    {
        Double diagonal=0.0;
        List_Process proceso=new List_Process();
        // a,b - the sequences, dw - the minimal distances matrix
        // w - the warping path
        //System.out.println(a.length+" "+b.length);
        int n=a.length,m=b.length,aux=0;
        
        Double d[][]=new Double[n][m]; // the euclidian distances matrix
        for(int i=0;i<n;i++)
        {
          for(int j=0;j<m;j++)
            {
//                if(i==j)
//                {
//                    System.out.println(distancia=distancia+proceso.distanciaHammit(a[i], b[j]));
//                }
                //d[i][j]=Math.abs(a[i]-b[j]);
                d[i][j]=proceso.distanciaHammit(a[i], b[j]);
            }
        }
        // determinate of minimal distance
        dw[0][0]=d[0][0];
        for(int i=1;i<n;i++)
        {
           dw[i][0]=d[i][0]+dw[i-1][0]; 
        }            
        for(int j=1;j<m;j++)
        {
            dw[0][j]=d[0][j]+dw[0][j-1];
        }            

        for(int i=1;i<n;i++) 
        {
            for(int j=1;j<m;j++) 
            {
                if(dw[i-1][j-1]<=dw[i-1][j]) 
                {
                    if(dw[i-1][j-1]<=dw[i][j-1])
                        if(i==j)//Priorizar diagonal
                        //System.out.println(dw[i][j]=d[i][j]+dw[i-1][j-1]-diagonal);
                            if(d[i][j]+dw[i-1][j-1]-diagonal<0)
                                dw[i][j]=-0.7;
                            else
                                dw[i][j]=d[i][j]+dw[i-1][j-1]-diagonal;                            
                        else
                            dw[i][j]=d[i][j]+dw[i-1][j-1];
                    else
                        if(i==j)//Priorizar diagonal
                        //System.out.println(dw[i][j]=d[i][j]+dw[i][j-1]-diagonal);
                            if(d[i][j]+dw[i][j-1]-diagonal<0)
                            dw[i][j]=-0.7;
                            else
                                dw[i][j]=d[i][j]+dw[i][j-1]-diagonal;                                
                        else
                            dw[i][j]=d[i][j]+dw[i][j-1];
                }
                else 
                {
                    if(dw[i-1][j]<=dw[i][j-1])
                        if(i==j)
                        //System.out.println(dw[i][j]=d[i][j]+dw[i-1][j]-diagonal);
                            if(d[i][j]+dw[i-1][j]-diagonal<0)
                                dw[i][j]=-1.0;
                            else
                                dw[i][j]=d[i][j]+dw[i-1][j]-diagonal;
                        else
                            dw[i][j]=d[i][j]+dw[i-1][j];
                    else
                        if(i==j)
                        //System.out.println(dw[i][j]=d[i][j]+dw[i][j-1]-diagonal);
                            if(d[i][j]+dw[i][j-1]-diagonal<0)
                                dw[i][j]=-1.0;
                            else
                                dw[i][j]=d[i][j]+dw[i][j-1]-diagonal;
                        else
                            dw[i][j]=d[i][j]+dw[i][j-1];
                }
            }
            if(i==400)
            {
                if(((100-(dw[i][m-1]/(double)m)*100))>75)
                {
                    /*System.out.println(100-(dw[i][m-1]/(double)m)*100+ "\t" +genero);
                    System.out.println("---i---" +i);*/
                    return 100-(dw[i][m-1]/(double)m)*100;
                }
                else
                {
                    return 0.0;
                }
            }
        }

        int i=n-1,j=m-1;
        double element=dw[i][j];
        //System.out.println("****Finalizado*******");
        //System.out.println(100-(element/(double)m)*100);
        return (100-(element/(double)m)*100);
        // determinate of warping path
        /*w.push(dw[i][j]);
        do{
            if(i>0&&j>0)
                if(dw[i-1][j-1]<=dw[i-1][j])
                    if(dw[i-1][j-1]<=dw[i][j-1])
                    {
                        i--;
                        j--;
                    } 
                    else 
                        j--;
                else
                    if(dw[i-1][j]<=dw[i][j-1])
                        i--; 
                    else j--;
            else 
                if(i==0)
                    j--; 
                else 
                    i--;
                    w.push(dw[i][j]);
        }
        while(i!=0||j!=0);
        while (!w.empty())
        {
            return (100-(element/(double)m)*100); 
        }
            
        //return element;
        return null;*/
    }
    
}
