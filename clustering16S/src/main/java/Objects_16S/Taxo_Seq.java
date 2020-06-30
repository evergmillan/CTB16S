/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Objects_16S;

/**
 *
 * @author cisei31
 */
public class Taxo_Seq implements Comparable<Taxo_Seq> {
    private String taxa;
    private String family;
    private String genre;
    private String id;
    private String sequence;
    private String protein_seq;
    private Double distance;
    private String grupo_genre;
    private int pos_genre;
    

    public Taxo_Seq() {
    }

    

    public String getTaxa() {
        return taxa;
    }

    public void setTaxa(String taxa) {
        this.taxa = taxa;
    }

    public String getGenre() {
        return genre;
    }

    public void setGenre(String genre) {
        this.genre = genre;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }
    
    public String getFamily() {
        return family;
    }

    public void setFamily(String family) {
        this.family = family;
    }
    
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getProtein_seq() {
        return protein_seq;
    }

    public void setProtein_seq(String protein_seq) {
        this.protein_seq = protein_seq;
    }
    
    public Taxo_Seq(String i) {
		this.sequence = i;
	}

    public Double getDistance() {
        return distance;
    }

    public void setDistance(Double distance) {
        this.distance = distance;
    }

    public String getGrupo_genre() {
        return grupo_genre;
    }

    public void setGrupo_genre(String grupo_genre) {
        this.grupo_genre = grupo_genre;
    }

    public int getPos_genre() {
        return pos_genre;
    }

    public void setPos_genre(int pos_genre) {
        this.pos_genre = pos_genre;
    }
    
    
    
	
    @Override
    public int compareTo(Taxo_Seq o) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
