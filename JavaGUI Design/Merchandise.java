package project;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.ResultSet;
import java.sql.Statement;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class Merchandise extends JFrame implements ActionListener{
	
	private JLabeledTextBox lb;
	private DataBase db;
	private ResultSet rs,rs1;
	private final int WINDOW_WIDTH = 550;
	private final int WINDOW_HEIGHT = 620;
	boolean modif=false;
	JListPanel l1;
	JListPanel l2,l3,l4;

	public Merchandise(DataBase db){
		super("Merchandise");
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
			rs = db.executeQuery("Select * from merchandise");
		 	//rs.first();
		 }
		catch(Exception ex){
			System.out.println(ex.getMessage());
		}	
		// Todo error:  Cannot read the array length because "Values" is null
		 moveFirst();		
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
		lb.addTextComponent("Reference", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Designation", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Unit price", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Type", JLabeledTextBox.TEXTFIELD);
		lb.setLabelsPreferredSize(80, 20);
		lb.setTextsPreferredSize(120, 20);
		lb.setBorder(BorderFactory.createTitledBorder("Information about the merchandise"));
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
		l1 = new JListPanel(5,20,100);
		l2 = new JListPanel(5,20,100);
		JButton select = new JButton(">>");
		select.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				l2.addItems(l1.getSelectedElements());
				l1.removeItems(l1.getSelectedElements());
			}
			});

		JButton retrieve = new JButton("<<");
		retrieve.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				l1.addItems(l2.getSelectedElements());
				l2.removeItems(l2.getSelectedElements());
			}
			});
			JPanel p = new JPanel(new GridLayout(2,1));
			p.add(select);
			p.add(retrieve);
			JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
			pan.add(l1);
			pan.add(p);
			pan.add(l2);
			
			pan.setBorder(BorderFactory.createTitledBorder("Raw material"));
		return pan;
	}
	private JPanel createList2(){
		l3 = new JListPanel(5,20,100);
		l4 = new JListPanel(5,20,100);
		JButton select = new JButton(">>");
		select.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				l4.addItems(l3.getSelectedElements());
				l3.removeItems(l3.getSelectedElements());
			}
			});
 
		JButton retrieve = new JButton("<<");
		retrieve.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				l3.addItems(l4.getSelectedElements());
				l4.removeItems(l4.getSelectedElements());
			}
			});
			JPanel p = new JPanel(new GridLayout(2,1));
			p.add(select);
			p.add(retrieve);
			JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
			pan.add(l3);
			pan.add(p);
			pan.add(l4);
			
			pan.setBorder(BorderFactory.createTitledBorder("Merchandise's Components"));
		return pan;
	}

	private void insert(){
		try{
			db.Insert("merchandise",lb.getValues());
			Object o[] = l4.getAllElements();
			for(int i=0;i<o.length;i++){
				db.executeUpdate("insert into compose_de values('"+lb.getValue(0)+"','"+o[i].toString()+"');");
			}
			Object o1[] = l2.getAllElements();
			for(int i=0;i<o1.length;i++){
				db.executeUpdate("insert into contient values('"+lb.getValue(0)+"','"+o1[i].toString()+"');");
			}
			JOptionPane.showMessageDialog(null,"Done","Production Management",JOptionPane.INFORMATION_MESSAGE);
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
		init();
	}
	private void delete(){
		try{
			Statement state =db.getConnection().createStatement();
			
			state.executeUpdate("delete from merchandise where reference='"+lb.getValue(0)+"'");
 
			JOptionPane.showMessageDialog(null,"Done","Production Management",JOptionPane.INFORMATION_MESSAGE);
			//moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void update(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("update merchandise set desig='"+lb.getValue(1)+"' where reference='"+lb.getValue(0)+"'");
			state.executeUpdate("update merchandise set unit_price='"+lb.getValue(2)+"' where reference="+lb.getValue(0));
			state.executeUpdate("update merchandise set type='"+lb.getValue(3)+"' where reference='"+lb.getValue(0)+"'");
			JOptionPane.showMessageDialog(null,"Done","Production Management",JOptionPane.INFORMATION_MESSAGE);
			///moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	
	private void updateLists() {
		try{
			l1.removeAllElements();
			Statement state =db.getConnection().createStatement();
			ResultSet res = state.executeQuery
					("SELECT code_mat from rawmaterial");
			//rs.first();
			while(res.next()){
				l1.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		try{
			l2.removeAllElements();
		}
		catch(Exception e){
			e.printStackTrace();
		}		
	
		try{
			l3.removeAllElements();
			Statement state =db.getConnection().createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,ResultSet.CONCUR_UPDATABLE);
			ResultSet res = state.executeQuery("SELECT code_mat FROM rawmaterial");
			//res.first();
			while(res.next()){
				l3.addItem(res.getString(1));
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		try{
			l4.removeAllElements();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Merchandise(db);
	}
}