package project;

import java.sql.*;
import java.util.Vector;

public class DataBase{
	
	public static final int ODBC = 1;
	public static final int MYSQL = 2;
	public static final int ORACLE=3;
	
	
	private  final String Accessdriver= "sun.jdbc.odbc.JdbcOdbcDriver";
	private  final String MYSQLdriver ="org.gjt.mm.mysql.Driver";
	private  final String ORACLEdriver="oracle.jdbc.driver.OracleDriver";	
	
	private  final String AccessBridge = "jdbc:odbc:";
	private  final String MysqlBridge = "jdbc:Mysql:";
	
	private Connection connect=null;
	private String driver;
	private String bridge;
	
	//////////////////////////////////////////////////:
	//Constructeur par defaut
	public DataBase(){
		driver = Accessdriver;
		bridge = AccessBridge;
	}
	////////////////////////////////////////////////::
	
	public void setDriver(int type){
		switch(type){
			case ODBC : driver = Accessdriver;
						bridge = AccessBridge;
						break;
			case MYSQL: driver = MYSQLdriver;
						bridge = MysqlBridge;
						break; 	
			case ORACLE:driver =ORACLEdriver;
						break;									
		}
	}
	
	//Constructeur avec le type du driver
	 public DataBase(int type){
	 	setDriver(type);
	 	
	 }
	 /////////////////////////////:::
	  public void Information(){
	 	
	 	System.out.println("Information sur la base de données : \n");
	 	try{
	 	
	 	DatabaseMetaData db = connect.getMetaData();
	 	System.out.println("__________||________________________");
	 	System.out.println("Connexion : "+ db.getURL());
	 	System.out.println("Driver    : "+ db.getDriverName());
	 	System.out.println("Version   : "+ db.getDriverVersion()); 	
	 	System.out.println("______________________________________");
	 	}
	 	catch(Exception e){
	 		System.out.println(e.getMessage());
	 	}
	 	
	 }
	 /////////////////////////////////////////////////////////////////
	 
	 public void prepare() {
	}
/////////////////////////////////////////////////////////////////////////////////////	 
	 
	 public void openDatabase(String DataBaseName){ 
	 	
	String url = bridge + DataBaseName;
	 	
	 	try{
	 		Class.forName(driver);
	 		connect = DriverManager.getConnection(url);
	 		System.out.println("base de donnees ouverte");
	 		prepare();
	 	}
	 	catch(Exception e){
	 		System.out.println("Erreur d'ouverture de la base de donnees \n"+e.getMessage());
	 	}
	 }
	 
/////////////////////////////////////////////////////////////////////////////////	 
	 
	 public void openDatabase(String DataBaseName,String pwd){
	 	String	url = bridge + DataBaseName;
	 	try{
	 		Class.forName(driver);
	 		connect = DriverManager.getConnection(url,"",pwd);
	 		prepare();
	 	}
	 	catch(Exception e){
	 		System.out.println(e.getMessage());
	 	}
	 }
	 
	 
	 
	 public void openDatabase(String DataBaseName,String login,String pwd){
	 	String url= bridge + DataBaseName;
	 	try{
	 		Class.forName(driver);
	 		connect = DriverManager.getConnection(url,login,pwd);
	 		prepare();
			}
	 	catch(Exception e){
	 		System.out.println(e.getMessage());
	 	}	
	 	
	 }
	 
	 public Connection getConnection(){
	 	return connect;
	 }
////////////////////////////requette SQL //////////////////////////////////////
//¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	 
	 public void Insert(String table, String[] values){
	 	String sql ="Insert into "+ table+ " values('";
		 	for(int i=0;i<values.length-1;i++){
		 		sql = sql + values[i] +"','";
		 	}
		 sql =sql + values[values.length-1]+"')";
		 
		 
		 try{
		 	Statement state = connect.createStatement();
		 	int r=state.executeUpdate(sql);
		 	state.close();
		 //	ResultSet resultat =state.executeUpdate(sql);
		 	
		 	
			 }
		catch(Exception e){
	 		System.out.println("Erreur d'insertion : " +e.getMessage());
	 		}	
	 	
	 }
	 
	 //-------------------------------------------------
	  public Vector select(String table){
	 	String sql="SELECT * FROM "+ table;
	 	Vector result = new Vector();
	 		try{
	 		Statement state = connect.createStatement();
	 		ResultSet rs= state.executeQuery(sql);
	 			String []T=null;

	 	while(rs.next()){
	 		int n= getColumnCount(rs);
	 	    T = new String[n];
	 		for(int i=0;i<n;i++)
	 			T[i]=rs.getString(i+1); 
	 		result.add(T);
	 	}
	 		return result;
	 		}
	 	catch(Exception e){
	 		System.out.println("Erreur de selection"+e.getMessage());
	 		return null;
	 		
	 		}	
	 		
	 	}
	 //---------------------------------------------------
	 	public int getColumnCount(ResultSet rs){
	 		
	 		try{ 
	 		ResultSetMetaData rsmd = rs.getMetaData();
	 	
	 		return rsmd.getColumnCount();
	 		}
	 		catch(Exception e){
	 		System.out.println(e.getMessage());
	 		return 0;	
	 	}
	 	}
	 //---------------------------------------------------
	 public int getColumnCount(String table){
	 	String sql="Select * from "+ table;
	 	try{ 
	 		Statement state = connect.createStatement();
	 		ResultSet resultat= state.executeQuery(sql);
	 		
	 		ResultSetMetaData rsmd = resultat.getMetaData();
	 		
	 		return rsmd.getColumnCount();
	 	}
	 	catch(Exception e){
	 		System.out.println(e.getMessage());
	 		return 0;
	 		}	
	 	
	 }
	 //------------------------------------------------------
	 
	 public String[] getColumnName(ResultSet rs){
	 	String T[]=null;
	 	
	 		try{
	 			ResultSetMetaData rsmd = rs.getMetaData();
	 			int n= rsmd.getColumnCount();
	 			T=new String[n]; 
	 			
	 			for(int i=0;i<n;i++){
	 				T[i]= rsmd.getColumnName(i+1);
	 				
	 			}
	 			return T;
	 		}	 			
	 	catch(Exception e){
	 		System.out.println(e.getMessage());
	 		return null;	
	 	}
	 	
	 }
	 
	 
	 
	 //------------------------------------------------------
	 public String[] getColumnName(String table){
	 	String sql="Select * from "+ table;
	 
	 	try{
	 		Statement state = connect.createStatement();
	 		ResultSet resultat= state.executeQuery(sql);
	 		
	 		return getColumnName(resultat);
	 	
	 	
	 		}
	 		
	 catch(Exception e){
	 	
	 		System.out.println(e.getMessage());
	 			return null;
	 		}
	 	//	return null;
	 }
	 //----------------------------------------------------
	 	public int getRowCount(ResultSet rs){
	 		int n=0;
	 		try{
	 			while(rs.next()){
	 				n++;
	 			}
	 		return n;	
	 		}
	 		catch(Exception e){
		 		System.out.println(e.getMessage());
		 			return 0;	
	 		}
	 	}
	 //----------------------------------------------------------
	 public int getRowCount(String table){
	 	String sql="Select * from "+ table;
	 	 		try{
	 		Statement state = connect.createStatement();
	 		ResultSet rs= state.executeQuery(sql);
	 			return getRowCount(rs);
	 		
	 	}
	 		catch(Exception e){
		 		System.out.println(e.getMessage());
		 			return 0;	
	 		}
	 	
	 	
	 	}
	 //----------------------------------------------------------
	  	public String [] getNextRecord(ResultSet rs) {
		try {
			String TC[] = null;
			rs.next();
			if(rs.isAfterLast()){
				return getLastRecord(rs);
			}
			else{
				int n = getColumnCount(rs);
				TC = new String[n];
				for (int i=0; i<n ; i++) {
				TC[i] = rs.getString(i+1);
				}
				return TC;
			}
		}
		catch(Exception e) {
			System.out.println("Erreur getNextRecord :"+e.getMessage());
			return null;
		}
	}
	  //----------------------------------------------------------  
	  	public String[] getNextRecord(String table){
		String sql="Select * from "+table;
		try{
			Statement state = connect.createStatement();
			ResultSet rs = state.executeQuery(sql);
			
			 	return getNextRecord(rs);
			 
		}
			catch(Exception e) {
			System.out.println( e.getMessage());
			return null;
		}
	}
	  
 //----------------------------------------------------------
	  	public String [] getPreviousRecord(ResultSet rs) {
		try {
			String TC[] = null;
			
			rs.previous();
			if(rs.isBeforeFirst()){
				return getFirstRecord(rs);
			}
			else{
				int n = getColumnCount(rs);
				TC = new String[n];
				for (int i=0; i<n; i++) {
				TC[i] = rs.getString(i+1);
				}
				return TC;	
			}
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
			return null;
		}
	}
	  //----------------------------------------------------------  
	  	public String[] getPreviousRecord(String table){
		String sql="Select * from "+table;
		try{
			Statement state = connect.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,
			ResultSet.CONCUR_UPDATABLE);
			ResultSet rs = state.executeQuery(sql);
			
			 	return getPreviousRecord(rs);
			 
		}
			catch(Exception e) {
			System.out.println( e.getMessage());
			return null;
		}
	}
	  

	  
	  
	  
	 //----------------------------------------------------------  
	  public String[] getLastRecord(ResultSet rs){
	  	String[] T =null;
	  	try{
	  		if(rs.last()){
	  			int n = getColumnCount(rs);
	  			 T=new String[n];
	  			for(int i=0;i<T.length;i++)
	  				T[i]= rs.getString(i+1);
	  		}
	  		return T;
	  	}
	  	catch(Exception e) {
			System.out.println("Erreur getNextRecord :"+e.getMessage());
			return null;
		}
	  	
	  }
	  //-----------------------------------------------------
	public String[] getLastRecord(String table){
		String sql="Select * from "+table;
		try{
			Statement state = connect.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,
			ResultSet.CONCUR_UPDATABLE);
			ResultSet rs = state.executeQuery(sql);
			
			 	return getLastRecord(rs);
			 
		}
			catch(Exception e) {
			System.out.println("Erreur getFirstRecord : " + e.getMessage());
			return null;
		}
	}
	  

	//---------------------------------------------------
	public String [] getFirstRecord(ResultSet rs) {
		try {
			String TC[]=null;
				
			if(rs.first()){
					
					int n = getColumnCount(rs);
					
					TC = new String[n];
					for (int i=0; i<n ; i++) {
						TC[i] = rs.getString(i+1);
					}
				}
			return TC;
		}
		catch(Exception e) {
			System.out.println("Erreur getFirstRecord : " + e.getMessage());
			return null;
		}
	}

	
	//---------------------------------------------------
	
	public String [] getFirstRecord(String table){
		String sql="Select * from "+table;
		try{
			Statement state = connect.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet rs = state.executeQuery(sql);
			
			 	return getFirstRecord(rs);
			 
		}
			catch(Exception e) {
			System.out.println("Erreur getFirstRecord : " + e.getMessage());
			return null;
		}
	}
	
	//----------------------------------------------------------
	
	 //----------------------------------------------------------
	
	public ResultSet executeQuery(String sql){
		
		try{
			Statement state = connect.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet rs = state.executeQuery(sql);
			return rs;
		}
			catch(Exception e) {
			System.out.println( e.getMessage());
			return null;
		}
	}
	public boolean executeUpdate(String sql){
		
		try{
			Statement state = connect.createStatement();
			state.executeUpdate(sql);
			return true;
		}
			catch(Exception e) {
			System.out.println( e.getMessage());
			return false;
		}
		
		
	}
	
 //----------------------------------------------------------	
	 		
}