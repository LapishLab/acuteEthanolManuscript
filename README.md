This code is to analyze the rat in vivo electrophysiology data for the acute ethanol manuscript

ratAcuteEthanol.m graphs each groups mean firing rates with error bars, 
and then correlates each neuron with each animal's respective intake, 
and graphs this data for Ensure only and Ensure + EtOH

ratRNNAcuteEthanol.m does exactly what ratAcuteEthanol.m does, 
except it doesnt generate a mean FR graphs, and correlates each 
neuron with each animal's respective brain EtOH predicted from 
the beiRNN instead of from their intake

beiRNNModel generates/trains the model fresh on input data of choice

beiRNNInstanceRuns runs the model with the model instances we already generated/tested
