package data;

import java.io.Serializable;

public class CascadeEvent implements Comparable<CascadeEvent>, Serializable{
    
    public int node;
    public double timestamp;

   
    public CascadeEvent(int nodeId,double timestamp) {
        super();
        this.node = nodeId;
        this.timestamp = timestamp;
    }

   

    public int getNodeId() {
        return node;
    }

    public void setNodeId(int nodeId) {
        this.node = nodeId;
    }

  
    public double getTimestamp(){
        return timestamp;
    }

   

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + node;
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
        CascadeEvent other = (CascadeEvent) obj;
        if (node != other.node)
            return false;
        return true;
    }

    @Override
    public int compareTo(CascadeEvent o) {
        if (timestamp < o.timestamp)
            return -1;
        if (timestamp > o.timestamp)
            return 1;
        return 0;
    }

}
