package data;

import java.io.Serializable;

public class CascadeEvent implements Comparable<CascadeEvent>, Serializable {

	private static final long serialVersionUID = -2498782211738252960L;

	public int node;
	public double timestamp;

	public CascadeEvent(int nodeId, double timestamp) {
		node = nodeId;
		this.timestamp = timestamp;
	}

	@Override
	public int compareTo(CascadeEvent o) {
		if (timestamp < o.timestamp)
			return -1;

		if (timestamp > o.timestamp)
			return 1;

		return 0;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;

		if (obj == null)
			return false;

		if (getClass() != obj.getClass())
			return false;

		final CascadeEvent other = (CascadeEvent) obj;

		if (node != other.node)
			return false;

		return true;
	}

	public int getNodeId() {
		return node;
	}

	public double getTimestamp() {
		return timestamp;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + node;
		return result;
	}

	public void setNodeId(int nodeId) {
		node = nodeId;
	}
}