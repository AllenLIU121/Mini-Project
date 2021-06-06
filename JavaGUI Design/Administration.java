package project;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.sql.*;

public class Administration extends JFrame implements ActionListener{
	
	private JLabeledTextBox lb;
	private DataBase db;
	private ResultSet rs;
	boolean modif=false;
	
	public Administration(DataBase db){
		super("Administration");
		lb = new JLabeledTextBox();	
		JLabel l = new JLabel("Administration");
		JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
		pan.add(l);
		getContentPane().add(pan,BorderLayout.NORTH);
		addButtons();
		initForm();
		try{
			this.db = db;
			init();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		setVisible(true);
		pack();
	}
	private void init(){
		 try{
		 	rs = db.executeQuery("Select * from users");
		 	//rs.first();
		 }
		catch(Exception ex){
			System.out.println(ex.getMessage());
	    }	
		 moveFirst();		
	}
	private void moveNext(){
		lb.setValues(db.getNextRecord(rs));
	}
	private void moveFirst(){
		lb.setValues(db.getFirstRecord(rs));
	}
	private void movePrevious(){
		lb.setValues(db.getPreviousRecord(rs));
	}
	private void moveLast(){
		lb.setValues(db.getLastRecord(rs));
	}
	private void addButtons(){
		String []s={"New","Save","Edit","Delete","Exit"};
		String []s1={"First","Previous","Next","Last"};
		JButtonImage bp1= new JButtonImage(s,1,s.length);
		JButtonImage bp2= new JButtonImage(s1,1,s1.length);
		bp1.addActionListener(this);
		bp2.addActionListener(this);
		
		JPanel p1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
		JPanel p2 = new JPanel(new FlowLayout(FlowLayout.CENTER));
		p1.add(bp1);
		p2.add(bp2);
		
		JPanel p3 = new JPanel(new GridLayout(2,1,2,2));
		p3.add(p1);
		p3.add(p2);
		
		getContentPane().add(p3,BorderLayout.SOUTH);
	}
	public void actionPerformed(ActionEvent ev){
		JButton b = (JButton)ev.getSource();
		if(b.getName().equals("Exit")){
			setVisible(false);
		}
		if(b.getName().equals("New")){
			lb.clear();
		}
		if(b.getName().equals("Save")){
			insert();
		}
		if(b.getName().equals("Delete")){
			delete();
		}
		if(b.getName().equals("Edit")){
			update();
		}
		if(b.getName().equals("First")){
			moveFirst();
		}
		if(b.getName().equals("Previous")){
			movePrevious();
		}
		if(b.getName().equals("Next")){
			moveNext();
		}
		if(b.getName().equals("Last")){
			moveLast();
		}
	}
	private void insert(){
		try{
			db.Insert("users",lb.getValues());
			JOptionPane.showMessageDialog(null,"Done","Production Management",JOptionPane.INFORMATION_MESSAGE);
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void delete(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("delete from users where id='"+lb.getValue(0)+"'");
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			moveNext();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void update(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("update users set pwd='"+lb.getValue(1)+"' where id='"+lb.getValue(0)+"'");
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			moveNext();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void initForm(){
		lb.addTextComponent("Identifier",JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Password",JLabeledTextBox.PASSWORD);
		lb.setTextsPreferredSize(120,20);
		//lb.setLockedAll(true);
		getContentPane().add(lb,BorderLayout.CENTER);
	}
	
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Administration(db);
	}
}