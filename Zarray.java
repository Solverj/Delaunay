
public class Zarray {

	int n;

	public void bygg(int[]x, int[]y, int[]z) {
		this.n = z.length;
		for (int i = 0; i < n; i++) {

			z[i] = (int) ((Math.cos(x[i]) * Math.sin(y[i])) * 100);
			System.out.println(z[i]);
		}
	}

}
