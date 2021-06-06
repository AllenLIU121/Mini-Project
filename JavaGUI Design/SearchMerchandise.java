package project;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableRowSorter;

import java.sql.*;

public class SearchMerchandise extends JFrame implements KeyListener{
	private JLabel labelSearch;
	private JPanel panel;
	private JPanel panelQuery;
	private JFlexGrid fg;
	private JTextField tf;
	DataBase db;

	public SearchMerchandise(DataBase db){
		super("Search Merchandise");
		this.db=db;
		tf= new JTextField(20);
		initJFlexGrid();
		initView();
		setVisible(true);
		pack();
	}
	
	private void initView(){
		panel = new JPanel();
		panelQuery = new JPanel(new FlowLayout(FlowLayout.LEFT));
		labelSearch = new JLabel("Search");

		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		fg.setBorder(BorderFactory.createTitledBorder("Result"));
		panelQuery.setBorder(BorderFactory.createTitledBorder("Query"));
		panelQuery.add(labelSearch);
		panelQuery.add(tf);
		tf.addKeyListener(this);
		panel.add(panelQuery);
		panel.add(fg);
		getContentPane().add(panel, BorderLayout.CENTER);
		setVisible(true);
		pack();
	}
	
	private void initJFlexGrid() {
		try{
			ResultSet rs = db.executeQuery("select * from merchandise;");
			int j=1;
			fg = new JFlexGrid(db.getRowCount(rs),db.getColumnCount(rs), true);
			fg.fill(rs);
			String s[]={"Reference","Designation","Unit Price","Type"};
			fg.setTitles(s);
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void keyTyped(KeyEvent e) {
		fg.search(tf.getText());
	}

	@Override
	public void keyPressed(KeyEvent e) {
		fg.search(tf.getText());
	}

	@Override
	public void keyReleased(KeyEvent e) {
		fg.search(tf.getText());
	}
	public static void main(String []args){
		DataBase db = new DataBase(DataBase.MYSQL);
		db.openDatabase("//localhost:3306/store management","root","");
		new SearchMerchandise(db);
	}
}