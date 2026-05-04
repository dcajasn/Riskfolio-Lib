###########################
Adanos Sentiment Functions
###########################

This module has optional helpers for using Adanos Market Sentiment API data
with Riskfolio-Lib portfolio workflows. The API covers stock sentiment from
Reddit, X / FinTwit, News and Polymarket, and can be used to build absolute
Black Litterman views.

The integration does not add a hard dependency. It uses Python's standard
library for the request and returns pandas objects that can be passed to the
existing Riskfolio-Lib estimators and portfolio classes.

Examples
========

::

    import riskfolio as rp

    assets = ["AAPL", "MSFT", "NVDA"]
    sentiment = rp.adanos_sentiment(
        assets,
        source="news",
        api_key="your-adanos-api-key",
    )
    P, Q = rp.adanos_sentiment_views(sentiment, assets, scale=0.02, threshold=0.1)

    # P and Q can be passed to Portfolio.blacklitterman_stats.

Module Functions
================

.. automodule:: SentimentFunctions
   :members:
   :private-members:
