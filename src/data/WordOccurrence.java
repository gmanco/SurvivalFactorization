package data;

import java.io.Serializable;

public class WordOccurrence implements Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    public int word;
    public int cnt;
    
   
    
    
    public WordOccurrence(int word, int cnt) {
        super();
        this.word = word;
        this.cnt = cnt;
    }

    
    



    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + word;
        return result;
    }






    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        WordOccurrence other = (WordOccurrence) obj;
        if (word != other.word)
            return false;
        return true;
    }






    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub

    }

}
