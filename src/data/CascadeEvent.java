package data;

import java.io.Serializable;

public class CascadeEvent implements Comparable<CascadeEvent>, Serializable{
    
    int nodeId;
    int cascadeId;
    double timestamp;

    public CascadeEvent(int nodeId, int cascadeId) {
        this.nodeId = nodeId;
        this.cascadeId = cascadeId;
        this.timestamp = -1;
    }

    public CascadeEvent(int nodeId, int cascadeId, double timestamp) {
        super();
        this.nodeId = nodeId;
        this.cascadeId = cascadeId;
        this.timestamp = timestamp;
    }

   

    public int getNodeId() {
        return nodeId;
    }

    public void setNodeId(int nodeId) {
        this.nodeId = nodeId;
    }

    public int getCascadeId() {
        return cascadeId;
    }

    public void setCascadeId(int cascadeId) {
        this.cascadeId = cascadeId;
    }
    
    public double getTimestamp(){
        return timestamp;
    }

   

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + cascadeId;
        result = prime * result + nodeId;
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
        if (cascadeId != other.cascadeId)
            return false;
        if (nodeId != other.nodeId)
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
