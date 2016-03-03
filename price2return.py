#!/usr/bin/env python

import pandas as pd
import numpy as np 

def compute_return(x, span):
	b = np.array(x[span: ])
	a = np.array(x[0 : len(x) - span])
	ratio = b / a
	return ratio


if __name__ == "__main__":

	df = pd.read_csv("pricedata.csv")

	ticker = list(df.columns.values)
	ticker.pop(0)

	nrow = len(df["Date"])
	ncol = len(ticker)
	span = 10 ## 10 business days

	data = np.empty([nrow - span, ncol])
	df_r = pd.DataFrame(data, columns = ticker)

	for column in df:
		if column != "Date":
			df_r[column] = compute_return(df[column], span)

	df_r.to_csv("return.csv", index = False, float_format="%.4f")
