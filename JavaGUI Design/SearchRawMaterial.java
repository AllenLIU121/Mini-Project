package project;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.sql.*;

public class SearchRawMaterial extends JFrame implements KeyListener{
	
	JFlexGrid fg;
	JTextField tf;
	DataBase db;
	public SearchRawMaterial(DataBase db){
		super("Search Raw Material");
		this.db=db;
		tf= new JTextField(20);
		initFlexGrid();
		
		JLabel l = new JLabel("Search:");
		JPanel p = new JPanel();

		p.setLayout(new BoxLayout(p,BoxLayout.Y_AXIS));
		JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		p1.add(l);
		p1.add(tf);
		tf.addKeyListener(this);
		
		p1.setBorder(BorderFactory.createTitledBorder("Query"));
		p.add(p1);
		fg.setBorder(BorderFactory.createTitledBorder("Results"));
		p.add(fg);
		getContentPane().add(p,BorderLayout.CENTER);
		setVisible(true);
		pack();
	}
	public void initFlexGrid(){
		try{
			ResultSet rs = db.executeQuery("select * from rawmaterial;");

			int j=1;
			fg = new JFlexGrid(db.getRowCount(rs),db.getColumnCount(rs), true);
			fg.fill(rs);
			String s[]={"Code","Designation"};
			fg.setTitles(s);
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		fg.search(tf.getText());
	}
	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		fg.search(tf.getText());
	}
	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		fg.search(tf.getText());
	}
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new SearchRawMaterial(db);
	}
	
}