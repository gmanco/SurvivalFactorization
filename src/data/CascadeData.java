package data;


import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.tint.impl.SparseIntMatrix2D;

public class CascadeData {
	   
        
        protected SparseDoubleMatrix2D cascadeEvents;
    
        /*
         * For each cascade it contains an ArrayList of cascade events (sorted)
         */
        protected ArrayList<CascadeEvent>[] cascadeEventsSrt; 
        
        /*
         * For each cascade it contains the set of active users
         */
        protected Set<Integer>[]activeNodesOnCascade;
        
        
        /*
         *  2-dim sparse matrix word x cascade
         */
        protected SparseIntMatrix2D CascadeContent; 
       
        /*
         * For each cascade it contains the set of active words
         */
        protected Set<Integer>[]activeWordsOnCascade;
         
        /*
         * Set of nodes ids
         */
        protected Set<Integer> nodes;
        
        /*
         * Set of cascadesIds
         */
        protected Set<Integer> cascades;
        
        /*
         * Set of word ids
         */
        protected Set<Integer> words;
        
        // properties
        protected int n_nodes;
        protected int n_cascades;
        protected int n_words = 0;
        protected double t_max;    
        
        
        public CascadeData(String file_events,String file_content){
            // here we read the file
            // file_cascade is in the format (u,c,t_u(c))
            // file_content is in the format (w,c,n_{w,c})
        		processEventFile(file_events);    
        	
            if (file_content != null)
            		processContentFile(file_content);            
        }
        
        
        
        /*
         * Return the activation timestamp for the pair
         * (node,cascade). Returns -1 if missing
         */
        public double getTimestamp(int nodeId,int cascadeId){
            Double t=cascadeEvents.get(nodeId, cascadeId);
            if(t==null){
                t=-1.0;
            }
            return t;
        }//getTimestamp
        
        
        public Set<Integer> activeNodesOnCascade(int cascadeId){
            return this.activeNodesOnCascade[cascadeId];
        }//activeNodesOnCascade
        
        
        public ArrayList<CascadeEvent> getSrtEventsForCascade(int cascadeId){
            return cascadeEventsSrt[cascadeId];
        }//getSrtEventsForCascade
        
        
        
        public Set<Integer> getInactiveNodesOnCascade(int cascadeId){
            HashSet<Integer> ris=new HashSet<Integer>(this.nodes);
            ris.removeAll(activeNodesOnCascade[cascadeId]);
            return ris;
        }//getInactiveNodesOnCascade
        
        
        public Set<Integer> activeWordsOnCascade(int cascadeId){
            return this.activeWordsOnCascade[cascadeId];
        }//activeNodesOnCascade
        
        
        public void getInfo(){
            System.out.format("Number of nodes: %d\n",n_nodes);
            System.out.format("Number of cascades %d\n",n_cascades);
            System.out.format("Number of words %d\n",n_words);
        }//getInfo        
                
        
        
        public double getNNodes(){
            return n_nodes;
        }
        
        public double getNCascades(){
            return n_cascades;
        }
        
        public double getNWords(){
            return n_words;
        }
                
        public double getTMax(){
            return t_max;
        }
        
       
        
        /*
         * Read Content file
         */
		private void processContentFile(String file_content) {
			try{
			    words=new TreeSet<Integer>();
			    
			    //read dimensions    
                BufferedReader br=new BufferedReader(new FileReader(file_content));
                String line=br.readLine();
                String tokens[];
                line=br.readLine();//skip header
                
                while(line!=null){
                    tokens=line.split("\t");
                    int word=Integer.parseInt(tokens[0]);
                    words.add(word);
                    
                }
                br.close();
                
                this.n_words=words.size();
                
                this.CascadeContent=new SparseIntMatrix2D(this.n_words,this.n_cascades);
                this.activeWordsOnCascade=new Set[this.n_cascades];
               
                br=new BufferedReader(new FileReader(file_content));
                line=br.readLine();
                line=br.readLine();//skip header
             
                Set<Integer> activeWordsTmp;
                while(line!=null){
                    tokens=line.split("\t");
                    int word=Integer.parseInt(tokens[0]);
                    int cascadeId=Integer.parseInt(tokens[1]);
                    int cnt=Integer.parseInt(tokens[2]);
                    this.CascadeContent.setQuick(word, cascadeId, cnt);
                    
                    activeWordsTmp=activeWordsOnCascade[cascadeId];
                    if(activeWordsTmp==null){
                        activeWordsTmp=new HashSet<Integer>();
                        activeWordsOnCascade[cascadeId]=activeWordsTmp;
                    }
                    activeWordsTmp.add(word);   
                }
                br.close();
			    
			}catch(Exception e){
			    e.printStackTrace();
			}
			
		}//processContentFile


		/*
		 * Read event file
		 */
		private void processEventFile(String file_events) {
		    try{
			
		        this.nodes=new TreeSet<Integer>();
	            this.cascades=new TreeSet<Integer>();
	            this.t_max=-1;
	            
    		    //read dimensions    
    		    BufferedReader br=new BufferedReader(new FileReader(file_events));
    			String line=br.readLine();
    			String tokens[];
    			line=br.readLine();//skip header
    			
    			while(line!=null){
    			    tokens=line.split("\t");
    			    int nodeId=Integer.parseInt(tokens[0]);
    			    int cascadeId=Integer.parseInt(tokens[1]);
    			    nodes.add(nodeId);
    			    cascades.add(cascadeId);
    			    line=br.readLine();
    			}
    			br.close();
			
    			//init variables
    			this.n_nodes=nodes.size();
    			this.n_cascades=cascades.size();
    			this.cascadeEvents=new SparseDoubleMatrix2D(n_nodes,n_cascades);
    			this.cascadeEventsSrt=new ArrayList[n_cascades];
    			this.activeNodesOnCascade=new HashSet[n_cascades];
    			
    			//read cascades
    			 br=new BufferedReader(new FileReader(file_events));
    	         line=br.readLine();
    	         line=br.readLine();//skip header
    	        
    	         ArrayList<CascadeEvent>tmp;
    	         Set<Integer>tmpSet;
    	         while(line!=null){
    	                tokens=line.split("\t");
    	                int nodeId=Integer.parseInt(tokens[0]);
    	                int cascadeId=Integer.parseInt(tokens[1]);
    	                double timestamp=Double.parseDouble(tokens[2]);
    	                CascadeEvent ce=new CascadeEvent(nodeId, cascadeId,timestamp);
    	              
    	                //set the timestamp
    	                cascadeEvents.set(nodeId, cascadeId, timestamp);
    	                
    	                // add event to CascadeEventsSrt
    	                tmp=cascadeEventsSrt[cascadeId];
    	                if(tmp==null){
    	                    tmp=new ArrayList<CascadeEvent>();
    	                    cascadeEventsSrt[cascadeId]=tmp;
    	                }
    	                tmp.add(ce);
    	              
    	                //add nodeId to activeUsersForCascade
    	                tmpSet=activeNodesOnCascade[cascadeId];
    	                if(tmpSet==null){
                            tmpSet=new HashSet<Integer>();
                            activeNodesOnCascade[cascadeId]=tmpSet;
                        }
    	                tmpSet.add(nodeId);
    	                
    	                //check t_max
    	                if(this.t_max<timestamp){
    	                    this.t_max=timestamp;
    	                }
    	                
    	          line=br.readLine();     
    	         }
    	         
    	         br.close();
    	         
    	         //sort events within cascades
    	         for(int cascadeId:cascades){
                     tmp=cascadeEventsSrt[cascadeId];
                     Collections.sort(tmp);        
                 }//sort 
    	         
    		    }catch(Exception e ){
    		        e.printStackTrace();
    		    }
		}//processEventFile

      
        

}
