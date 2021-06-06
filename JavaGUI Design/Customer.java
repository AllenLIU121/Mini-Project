package project;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.ResultSet;
import java.sql.Statement;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

public class Customer extends JFrame implements ActionListener{
	
	private JLabeledTextBox lb;
	private DataBase db;
	private ResultSet rs;
	private JLabel label;
	boolean edit=false;
	
	public Customer(DataBase db){
		super("Customer");
		initView();
		try{
			this.db = db;
			init();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		setVisible(true);
		pack();
		setResizable(false);
	}

	private void initView(){
		lb = new JLabeledTextBox();	
		addLabels();
		initForm();
		addButtons();
	}
	
	private void init(){
		 try{
			rs = db.executeQuery("Select * from customer");
			System.out.println("RS2:" + rs);
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

	private void addLabels(){
		label = new JLabel("Customer");
		label.setHorizontalAlignment(SwingConstants.CENTER);
		getContentPane().add(label, BorderLayout.NORTH);
	}

	private void addButtons(){
		String []s1={"New","Save","Edit","Delete","Exit"};
		String []s2={"First","Previous","Next","Last"};
		JButtonImage bp1= new JButtonImage(s1,1,s1.length);
		JButtonImage bp2= new JButtonImage(s2,1,s2.length);
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
	
	private void initForm(){
		lb.addTextComponent("Customer id", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Full name", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Address", JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Telephone", JLabeledTextBox.TEXTFIELD);
		lb.setLabelsPreferredSize(80, 20);
		lb.setTextsPreferredSize(120, 20);
		add(lb);
	}

	private void insert(){
		if(edit == false) {
			try{
				db.Insert("Customer",lb.getValues());
				JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			}
			catch(Exception e){
				JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
			}
			init();
		}
		else
			JOptionPane.showMessageDialog(null,"Duplicate Entry","Production Management",JOptionPane.ERROR_MESSAGE);

	}
	private void delete(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("delete from customer where customer_id ='"+lb.getValue(0)+"'");
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
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("update Customer set full_name_cs='"+lb.getValue(1)+"' where customer_id ='"+lb.getValue(0)+"'");
			state.executeUpdate("update Customer set address_cs='"+lb.getValue(2)+"' where customer_id ='"+lb.getValue(0)+"'");
			state.executeUpdate("update Customer set telephone_cs='"+lb.getValue(3)+"' where customer_id ='"+lb.getValue(0)+"'");
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			///moveNext();
			init();
			edit = false;
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}

	
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Customer(db);
	}
}