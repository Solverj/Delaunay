class shuftTest{

    public static void main(String[] args){
	int xshift = 3;
	int maxXY = 27;
        while ((1<<xshift) < maxXY) xshift++;
	System.out.println(xshift);
	System.out.println(1<<4);
    }

}
