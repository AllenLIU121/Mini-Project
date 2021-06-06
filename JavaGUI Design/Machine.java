package project;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.sql.*;

public class Machine extends JFrame implements ActionListener{
	
	private JLabeledTextBox lb;
	private DataBase db;
	private ResultSet rs;
	private final int WINDOW_WIDTH = 550;
	private final int WINDOW_HEIGHT = 620;
	boolean modif=false;
	JListPanel list_merchandise1, list_merchandise2;
	JListPanel list_raw_material1, list_raw_material2;
	public Machine(DataBase db){
		super("Machine");
		initView();
		try{
			this.db = db;
		}
		catch(Exception e){
			e.printStackTrace();
		}
		init();
		
	}

	private void initView(){
		setPreferredSize(new Dimension(WINDOW_WIDTH, WINDOW_HEIGHT));
		setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
		lb = new JLabeledTextBox();	
		initForm();
		add(createList1());
		add(createList2());
		addButtons();	
		setVisible(true);
		pack();		
	}

	private void init(){
		 try{
			rs = db.executeQuery("Select * from machine");
		 	//rs.first();
		 }
		catch(Exception ex){
			System.out.println(ex.getMessage());
		}
		 moveFirst();		
		 updateLists();
	}
	private void moveNext(){
		lb.setValues(db.getNextRecord(rs));
		updateLists();
	}
	private void moveFirst(){
		lb.setValues(db.getFirstRecord(rs));
		updateLists();
	}
	private void movePrevious(){
		lb.setValues(db.getPreviousRecord(rs));
		updateLists();
	}
	private void moveLast(){
		lb.setValues(db.getLastRecord(rs));
		updateLists();
	}

	private void addButtons(){
		String []s1={"New","Save","Edit","Delete","Exit"};
		String []s2={"First","Previous","Next","Last"};
		JButtonImage bp1= new JButtonImage(s1,1,s1.length);
		JButtonImage bp2= new JButtonImage(s2,1,s2.length);
		bp1.addActionListener(this);
		bp2.addActionListener(this);
		JPanel p1 = new JPanel();
		JPanel p2 = new JPanel();
		p1.add(bp1);
		p2.add(bp2);
		JPanel p3 = new JPanel(new GridLayout(0,1,2,2));
		p3.add(p1);
		p3.add(p2);
		getContentPane().add(p3);
	}
	private void initForm(){
		lb.addTextComponent("Machine id", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Designation", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Description", JLabeledTextBox.TEXTFIELD);
		lb.setLabelsPreferredSize(80, 20);
		lb.setTextsPreferredSize(120, 20);
		lb.setBorder(BorderFactory.createTitledBorder("Information about the machine"));
		add(lb);
	}
	public void actionPerformed(ActionEvent ev){
		JButton button = (JButton)ev.getSource();
		switch(button.getName()){
			case "New":
				lb.clear();
				lb.setLockedAll(false);
				break;
			case "Save":
				insert();
				break;
			case "Edit":
				update();
				break;
			case "Delete":
				delete();
				break;
			case "Exit":
				setVisible(false);
				break;
			
			case "First":
				moveFirst();
				break;
			case "Previous":
				movePrevious();
				break;
			case "Next":
				moveNext();
				break;
			case "Last":
				moveLast();
			default:
				break;
		}
	}
	private JPanel createList1(){
		list_merchandise1 = new JListPanel(5,20,100);
		list_merchandise2 = new JListPanel(5,20,100);
		JButton select = new JButton(">>");
		select.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				list_merchandise2.addItems(list_merchandise1.getSelectedElements());
				list_merchandise1.removeItems(list_merchandise1.getSelectedElements());
			}
			});
		JButton retrieve = new JButton("<<");
		retrieve.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				list_merchandise1.addItems(list_merchandise2.getSelectedElements());
				list_merchandise2.removeItems(list_merchandise2.getSelectedElements());
			}
			});
			JPanel p = new JPanel(new GridLayout(2,1));
			p.add(select);
			p.add(retrieve);
			JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
			pan.add(list_merchandise1);
			pan.add(p);
			pan.add(list_merchandise2);
			
			pan.setBorder(BorderFactory.createTitledBorder("Raw Materials to Process"));
		return pan;
	}
	private JPanel createList2(){
		list_raw_material1 = new JListPanel(5,20,100);
		list_raw_material2 = new JListPanel(5,20,100);
		JButton select = new JButton(">>");
		select.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				list_raw_material2.addItems(list_raw_material1.getSelectedElements());
				list_raw_material1.removeItems(list_raw_material1.getSelectedElements());
			}
			});

		JButton retrieve = new JButton("<<");
		retrieve.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				list_raw_material1.addItems(list_raw_material2.getSelectedElements());
				list_raw_material2.removeItems(list_raw_material2.getSelectedElements());
			}
		});
		JPanel p = new JPanel(new GridLayout(2,1));
		p.add(select);
		p.add(retrieve);
		JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
		pan.add(list_raw_material1);
		pan.add(p);
		pan.add(list_raw_material2);
			
		pan.setBorder(BorderFactory.createTitledBorder("Merchandise"));
		return pan;
	}

	private void insert(){
		try{
			db.Insert("machine",lb.getValues());
			Object o[] = list_raw_material2.getAllElements();
			for(int i=0;i<o.length;i++){
				db.executeUpdate("insert into passed_by values('"+lb.getValue(0)+"','"+o[i].toString()+"');");
			}
			Object o1[] = list_merchandise2.getAllElements();
			for(int i=0;i<o1.length;i++){
				db.executeUpdate("insert into processed_by values('"+o1[i].toString()+"','"+lb.getValue(0)+"');");
			}
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
		init();
	}
	private void delete(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("delete from machine where code_m='"+lb.getValue(0)+"'");
			state.executeUpdate("delete from processed_by where code_m='"+lb.getValue(0)+"'");
			state.executeUpdate("delete from passed_by where code_m='"+lb.getValue(0)+"'");
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			//moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void update(){
		try{
			Statement state = db.getConnection().createStatement();
			state.executeUpdate("update machine set designation='"+lb.getValue(1)+"' where code_m='"+lb.getValue(0)+"'");
			state.executeUpdate("update machine set description='"+lb.getValue(2)+"' where code_m='"+lb.getValue(0)+"'");
			state.executeUpdate("delete from processed_by where code_m='"+lb.getValue(0)+"'");
			state.executeUpdate("delete from passed_by where code_m='"+lb.getValue(0)+"'");
	
			Object o[] = list_raw_material2.getAllElements();
			for(int i=0;i<o.length;i++){
				state.executeUpdate("insert into passed_by values('"+lb.getValue(0)+"','"+o[i].toString()+"');");
			}
			Object o1[] = list_merchandise2.getAllElements();
			for(int i=0;i<o1.length;i++){
				state.executeUpdate("insert into processed_by values('"+o1[i].toString()+"','"+lb.getValue(0)+"');");
			}
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			///moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}

	
	private void updateLists() {
		try{
			list_merchandise1.removeAllElements();
			Statement state =db.getConnection().createStatement();
			ResultSet res = state.executeQuery
					("SELECT distinct(passed_by.reference) from machine, passed_by WHERE machine.code_m=passed_by.code_m and  passed_by.code_m<>'"+lb.getValue(0)+"'");
			//rs.first();
			while(res.next()){
				list_merchandise1.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		try{
			list_merchandise2.removeAllElements();
			Statement state =db.getConnection().createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet res = state.executeQuery
					("SELECT distinct(reference) from passed_by WHERE code_m='"+lb.getValue(0)+"'");
			//res.first();
			while(res.next()){
				list_merchandise2.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}		
	
		try{
			list_raw_material1.removeAllElements();
			Statement state =db.getConnection().createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet res = state.executeQuery("SELECT distinct(machine.code_m) FROM machine, processed_by "
			+ "WHERE machine.code_m=processed_by.code_m and processed_by.code_m<>'"+lb.getValue(0)+"'");
			//res.first();
			while(res.next()){
				list_raw_material1.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		try{
			
			list_raw_material2.removeAllElements();
			Statement state =db.getConnection().createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet res = state.executeQuery("SELECT distinct(code_mat) FROM processed_by WHERE " + 
					"code_m='"+lb.getValue(0)+"'");
			//res.first();
			while(res.next()){
				list_raw_material2.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Machine(db);
	}
}