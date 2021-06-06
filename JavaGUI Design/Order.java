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
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class Order extends JFrame implements ActionListener{
	
	private JLabeledTextBox lb;
	private DataBase db;
	private ResultSet rs;
	JListPanel l1,l2;
	JComboBox b;
	public Order(DataBase db){
		super("Order");
	
		lb = new JLabeledTextBox();	
		JLabel l = new JLabel("Order");
		b =  new JComboBox();
		JPanel pan = new JPanel(new FlowLayout(FlowLayout.CENTER));
		pan.add(l);
		getContentPane().add(pan,BorderLayout.NORTH);
		addButtons();
		try{
			this.db = db;
		}
		catch(Exception e){
			e.printStackTrace();
		}
		initForm();
		init();
		setVisible(true);
		pack();
	}
	private void init(){
		 try{
		 	rs = db.executeQuery("Select * from invoiced_order");
		 	//rs.first();
		 }
		catch(Exception ex){
			System.out.println(ex.getMessage());
	    }	
		 moveFirst();		
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
		
		JPanel p3 = new JPanel(new GridLayout(0,1,2,2));
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
			lb.setLockedAll(false);
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
	private void moveNext(){
		String []s= db.getNextRecord(rs);
		b.setSelectedItem(s[1]);
		String s1[]={s[0],s[2],s[3],s[4],s[5]};
		lb.setValues(s1);
	}
	private void moveFirst(){
		String []s= db.getFirstRecord(rs);
		b.setSelectedItem(s[1]);
		String s1[]={s[0],s[2],s[3],s[4],s[5]};
		lb.setValues(s1);
	}
	private void movePrevious(){
		String []s= db.getPreviousRecord(rs);
		b.setSelectedItem(s[1]);
		String s1[]={s[0],s[2],s[3],s[4],s[5]};
		lb.setValues(s1);
	}
	private void moveLast(){
		String []s= db.getLastRecord(rs);
		b.setSelectedItem(s[1]);
		String s1[]={s[0],s[2],s[3],s[4],s[5]};
		lb.setValues(s1);
	}
	private void insert(){
		try{
			String []s={lb.getValue(0),b.getSelectedItem().toString(),lb.getValue(1),
							lb.getValue(2),lb.getValue(3),lb.getValue(4)};
			
			//init();
			db.Insert("invoiced_order",s);
			Object o[]=l2.getAllElements();
			for(int i=0;i<o.length;i++){
				db.executeUpdate("insert into compose values('"+lb.getValue(0)+"','"+o[i].toString()+"')");
			}
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			init();
		}
		catch(Exception e){
		}
	}
	private void delete(){
		try{
			Statement state =db.getConnection().createStatement();
			state.executeUpdate("delete from invoiced_order where ncde="+lb.getValue(0));
			state.executeUpdate("delete from composer where ncde="+lb.getValue(0));
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			//moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private JPanel createList(){
		l1 = new JListPanel(5,20,100);
		l2 = new JListPanel(5,20,100);
		JButton select = new JButton(">>");
		select.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent ev){
				l2.addItems(l1.getSelectedElements());
				l1.removeItems(l1.getSelectedElements());
			}
			});
			fill1();
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
			
			pan.setBorder(BorderFactory.createTitledBorder("Merchandise"));
		return pan;
	}

	private void update(){
		try{
			db.executeUpdate("update invoiced_order set date_of_order='"+lb.getValue(1)+"' where ncde="+lb.getValue(0));
			db.executeUpdate("update invoiced_order set date_of_delivery='"+lb.getValue(2)+"' where ncde="+lb.getValue(0));
			db.executeUpdate("update invoiced_order set payment_method='"+lb.getValue(3)+"' where ncde="+lb.getValue(0));
			db.executeUpdate("update invoiced_order set state='"+lb.getValue(4)+"' where ncde="+lb.getValue(0));
			JOptionPane.showMessageDialog(null,"Done successfully","Production Management",JOptionPane.INFORMATION_MESSAGE);
			///moveNext();
			init();
		}
		catch(Exception e){
			JOptionPane.showMessageDialog(null,e.getMessage(),"Production Management",JOptionPane.ERROR_MESSAGE);
		}
	}
	private void fill(){
		try{
			ResultSet r=db.executeQuery("select * from customer");
			r.first();
			b.addItem(r.getString(1));

			while(r.next()){
				b.addItem(r.getString(1));
			}
		}
		
		catch(Exception e){
			System.out.println(e.getMessage());
		}
	}
		private void fill1(){
		try{
			ResultSet r=db.executeQuery("select * from merchandise");
			r.first();
		l1.addItem(r.getString(1));

			while(r.next()){
				l1.addItem(r.getString(1));
			}
		}
		
		catch(Exception e){
			System.out.println(e.getMessage());
		}
	}

	private void initForm(){
		fill();
		JLabel l = new JLabel("Customer ID");
		JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		l.setPreferredSize(new Dimension(120,20));
		p1.add(l);
		p1.add(b);
		lb.add(p1);
		lb.addTextComponent("Order number",JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Date of order",JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Date of delivery",JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("Payement method",JLabeledTextBox.TEXTFIELD);
		lb.addTextComponent("State",JLabeledTextBox.TEXTFIELD);
		lb.setTextsPreferredSize(120,20);
		//lb.setLockedAll(true);
		lb.setLocked(0,true);
		lb.setLabelsPreferredSize(120,20);
		JPanel p = new JPanel();
		p.setLayout(new BoxLayout(p,BoxLayout.Y_AXIS));
		p.add(lb);
		p.add(createList());
		getContentPane().add(p,BorderLayout.CENTER);

	}
		public boolean exist(String req){
		if(db.getRowCount(req)>0){
			return true;
		}
		else
			return false;
	}
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Order(db);
	}
}