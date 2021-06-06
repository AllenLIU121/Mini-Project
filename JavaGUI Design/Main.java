package project;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JMenuItem;

public class Main extends JFrame implements ActionListener{
	
	Container content;
	Image img = Toolkit.getDefaultToolkit().getImage("Images/Production2.png");
	private DataBase db;
	public Main(DataBase db){
		super("Production management");
		this.db = db;
		content = getContentPane();	
		addJMenu();
		addButtons();
		DrawImage des = new DrawImage(img);
		des.repaint();
		content.add(des,BorderLayout.CENTER);
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		pack();
		setSize(500, 500);
		setVisible(true);
	}
	
	public void addJMenu(){
		String [][]s={{"File","Administration","Exit"},
					{"Consultation","Merchandise","Raw material","Customer","Machine","Store","Order"},
					{"Search","Search Merchandise","Search Raw material","Search Customer"}
					};
		SimpleMenu sp = new SimpleMenu(s, new MenuBarListener());
		setJMenuBar(sp);			
	}
	public void addButtons(){
		String s[]={"Raw material","Merchandise","Machine","Store","Customer","Order","Exit"};
		JToolBarButton bp = new JToolBarButton(s);
		bp.addActionListener(this);
		content.add(bp,BorderLayout.NORTH);
	}
	private class MenuBarListener implements ActionListener{
		public void actionPerformed(ActionEvent ev){
			JMenuItem it =(JMenuItem)ev.getSource();
			if(it.getText().equals("Exit")){
				System.exit(-1);
			}
			if(it.getText().equals("Administration")){
				new Administration(db).setLocation(140,45);
			}
			if(it.getText().equals("Raw material")){
				new RawMaterial(db).setLocation(140,45);
			}
			if(it.getText().equals("Machine")){
				new Machine(db).setLocation(140,45);
			}
			if(it.getText().equals("Store")){
				new Store(db).setLocation(140,45);
			}
			if(it.getText().equals("Customer")){
				new Customer(db).setLocation(140,45);
			}
			if(it.getText().equals("Merchandise")){
				new Merchandise(db).setLocation(140,45);
			}
			if(it.getText().equals("Order")){
				new Order(db).setLocation(140,45);
			}
			if(it.getText().equals("Search Merchandise")){
				new SearchMerchandise(db).setLocation(140,45);
			}
			if(it.getText().equals("Search Customer")){
				new SearchCustomer(db).setLocation(140,45);
			}
			if(it.getText().equals("Search Raw material")){
				new SearchRawMaterial(db).setLocation(140,45);
			}
		}
	}
	public void actionPerformed(ActionEvent ev){
		JButton b = (JButton)ev.getSource();
		if(b.getName().equals("Machine")){
			new Machine(db).setLocation(140,45);
		}
		if(b.getName().equals("Raw material")){
			new RawMaterial(db).setLocation(140,45);
		}
		if(b.getName().equals("Store")){
			new Store(db).setLocation(140,45);
		}
		if(b.getName().equals("Customer")){
			new Customer(db).setLocation(140,45);
		}
		if(b.getName().equals("Merchandise")){
			new Merchandise(db).setLocation(140,45);
		}
		if(b.getName().equals("Order")){
			new Order(db).setLocation(140,45);
		}
		if(b.getName().equals("Exit")){
			System.exit(-1);
		}
	}
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new Main(db);
	}
}