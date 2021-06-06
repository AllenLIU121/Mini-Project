package project;

import java.awt.*;
import java.sql.ResultSet;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableRowSorter;


public class JFlexGrid extends JPanel{
	
	private JTable table;
	private DefaultTableModel  tm;
	private TableRowSorter<DefaultTableModel> rowSorter;

	public JFlexGrid(int Row,int Col){
		super(new BorderLayout());
		
		tm = new DefaultTableModel(Row,Col);
		table = new JTable(tm);
		add("Center",new JScrollPane(table));
	}
	
	public JFlexGrid(int Row,int Col, boolean setSorter){
		super(new BorderLayout());
		
		tm = new DefaultTableModel(Row,Col);
		table = new JTable(tm);
		if(setSorter) {
			rowSorter = new TableRowSorter<DefaultTableModel>(tm);
			table.setRowSorter(rowSorter);
		}
		add("Center",new JScrollPane(table));
	}
	
	public void search(String str) {
		if(rowSorter!=null) {
			if (str.length() == 0) {
				rowSorter.setRowFilter(null);
            } 
			else {
				rowSorter.setRowFilter(RowFilter.regexFilter(str));
            }
        }
	}
	public void setHeaderTitle(int Col,String Title){
		table.getColumnModel().getColumn(Col).setHeaderValue(Title);
	}
	
	public void setTitles(String[] titles){
		for(int i=0;i<titles.length;i++){
			setHeaderTitle(i,titles[i]);
		}
	}
	
	public void fill(ResultSet resultSet) {
		try {
			int cols = resultSet.getMetaData().getColumnCount();
			resultSet.first();
			for(int i=0;i<cols;i++){
				setTextMatrix(0,i,resultSet.getString(i+1));
			}
			int j=1;
			while(resultSet.next()){
				for(int i=0;i<cols;i++){
					setTextMatrix(j,i,resultSet.getString(i+1));
				}
				j++;
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	public void setRowSorter(RowSorter rowSorter) {
		table.setRowSorter(rowSorter);
	}
	public void setTextMatrix(int Row,int Col,Object text){
		tm.setValueAt(text,Row,Col);		
	}
	
	public Object getTextMatrix(int Row,int Col){
		return tm.getValueAt(Row,Col);
	}
	
	public int getSelectCol(){
		return table.getSelectedColumn();
	}
	public int getSelectedRow(){
		return table.getSelectedRow();
	}
	
	
	//いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
	
	public static void main(String[] args){
		
		JFlexGrid fg = new JFlexGrid(2,3);
		fg.setTextMatrix(1,2,"OK");
		String s[]={"Nom","prenom","Adresse"};
		fg.setTitles(s);
		
		JFrame f =new JFrame("Table");
		f.add(fg);
		f.pack();
		f.setVisible(true);
	}
} 