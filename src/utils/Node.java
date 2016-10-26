package utils;

public class Node implements Comparable<Node> {
	private int id;
	private double value;

	public Node() {
		// do nothing!
	}

	public Node(int id, double value) {
		this.id = id;
		this.value = value;
	}

	@Override
	public int compareTo(Node o) {
		if (value < o.value)
			return -1;

		if (value > o.value)
			return 1;

		if (id < o.id)
			return -1;

		if (id > o.id)
			return 1;

		return 0;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Node)
			return id == ((Node) obj).id && value == ((Node) obj).value;

		return false;
	}

	public int getId() {
		return id;
	}

	public double getValue() {
		return value;
	}

	public void setId(int id) {
		this.id = id;
	}

	public void setValue(double value) {
		this.value = value;
	}

	@Override
	public String toString() {
		return "id: " + id + ", value: " + value;
	}
}