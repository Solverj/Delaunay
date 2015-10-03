public class beregnKote{

  int x1 = 0,
      x2 = 5,
      y1 = 1,
      y2 = 2;

  public beregnKote {

  }
  public static void main(String[] args) {
    new beregnKote();
    }
  double beregnDistanseMellomToPunkt(int x1, int y1, int x2, int y2){
    return (double) Math.sqrt((Math.abs((x2 - x1))^2) + (Math.abs((y2 - y1))^2));
  }
}
