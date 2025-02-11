######################################################
Riskfolio-XL: Riskfolio-Lib add-in for Microsoft Excel
######################################################

.. raw:: html

    <a href="https://www.kqzyfj.com/click-101359873-15150084?url=https%3A%2F%2Flink.springer.com%2Fbook%2F9783031843037" target="_blank">
        <button style="padding:10px 20px; font-size:16px; background-color: #FFA500; color:white; border:none; border-radius:5px; cursor:pointer;">
            Buy Advanced Portfolio Optimization Book on Springer
        </button>
    </a>
    <br>
    <br>

.. image:: https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86
 :target: https://github.com/sponsors/dcajasn

.. raw:: html
   
    <br>
   
.. raw:: html

    <a href='https://ko-fi.com/B0B833SXD' target='_blank'><img height='36'style='border:0px;height:36px;' src='https://cdn.ko-fi.com/cdn/kofi1.png?v=2' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>


Description
===========

Riskfolio-XL is a Microsoft Excel add-in based on `PyXLL <https://www.pyxll.com/index.html>`_ library, that allows users use the same features of Riskfolio-Lib in Excel through Riskfolio-XL spreadsheet functions. Its objective is to help non-programming users to build investment portfolios based on mathematically complex models with low effort and to support the maintenance and further development of Riskfolio-Lib.

Installation
============

Riskfolio-XL is only available on Windows and it requires a valid installation of PyXLL package and PyXLL add-in. To install PyXLL and PyXLL add-in, you can find the PyXLL installation instructions in the following `link <https://www.pyxll.com/docs/userguide/installation/firsttime.html>`_.

After installing the PyXLL package and PyXLL add-in, the latest stable release of Riskfolio-XL (and older versions) can be installed from PyPI:

  ::
    
      pip install riskfolio-xl


After installing the Riskfolio-XL package you will have access to the **TRIAL COPY** of Riskfolio-XL, this version is limited to work only with portfolios of 7 assets and risk factor models of 3 risk factors.

To access the **PURCHASED COPY** of Riskfolio-XL, you need a valid license. To get a Riskfolio-XL license you have to purchase it paying a monthly or annual subscription:

.. raw:: html

   <style>
   .tab {
   overflow: hidden;
   border: 1px solid #ccc;
   background-color: #f1f1f1;
   }

   /* Style the buttons inside the tab */
   .tab button {
   background-color: inherit;
   float: left;
   border: none;
   outline: none;
   cursor: pointer;
   padding: 14px 16px;
   transition: 0.3s;
   font-size: 17px;
   }

   /* Change background color of buttons on hover */
   .tab button:hover {
   background-color: #ddd;
   }

   /* Create an active/current tablink class */
   .tab button.active {
   background-color: #ccc;
   }

   /* Style the tab content */
   .tabcontent {
   display: none;
   padding: 6px 12px;
   border: 1px solid #ccc;
   border-top: none;
   }
   </style>
   <div class="tab">
   <button class="tablinks" onclick="openTab(event, 'paypal-container-DLABCNSPX8LZL')" id="defaultOpen">Monthly License</button>
   <button class="tablinks" onclick="openTab(event, 'paypal-container-KKZBWK6JK8ZDA')">Annual License</button>
   </div>

   <script src="https://www.paypal.com/sdk/js?client-id=BAA_FQBdhZjxYgI2N5DACAiN0--Lkv3sO9Kj0LKlFq9BWpNha13pFGIjK3X9qumuLmkh9oOPFdoSb-mJvc&components=hosted-buttons&disable-funding=venmo&currency=USD"></script>
   <div id="paypal-container-DLABCNSPX8LZL" class="tabcontent"></div>
   <div id="paypal-container-KKZBWK6JK8ZDA" class="tabcontent"></div>
   <script>
   paypal.HostedButtons({
      hostedButtonId: "DLABCNSPX8LZL",
   }).render("#paypal-container-DLABCNSPX8LZL")

   paypal.HostedButtons({
      hostedButtonId: "KKZBWK6JK8ZDA",
   }).render("#paypal-container-KKZBWK6JK8ZDA")

   function openTab(evt, tabId) {
   var i, tabcontent, tablinks;
   tabcontent = document.getElementsByClassName("tabcontent");
   for (i = 0; i < tabcontent.length; i++) {
      tabcontent[i].style.display = "none";
   }
   tablinks = document.getElementsByClassName("tablinks");
   for (i = 0; i < tablinks.length; i++) {
      tablinks[i].className = tablinks[i].className.replace(" active", "");
   }
   document.getElementById(tabId).style.display = "block";
   evt.currentTarget.className += " active";
   }

   document.getElementById("defaultOpen").click();
   </script>
   </br>

**After paying, you need to send us an email to** `orenji.eirl@gmail.com <orenji.eirl@gmail.com>`_ **and you will receive your Riskfolio-XL license (within 24 hours) and a discount code to purchase the PyXLL package.**

Then, you have to add your Riskfolio-XL license to the pyxll.cfg file:

* First, click on About PyXLL button of Riskfolio-XL add-in as shown in the image below:

.. image:: images/Images-XL/Installation_1.png

* Then, click on the config file link as shown in the image below:

.. image:: images/Images-XL/Installation_2.png

* Finally, write your Riskfolio-XL license in the riskfolio_xl_key parameter as shown in the image below:

.. image:: images/Images-XL/Installation_3.png

Citing
======

If you use Riskfolio-Lib for published work, please use the following BibTeX entry:

::

    @misc{riskfolioxl,
          author = {Dany Cajas},
          title = {Riskfolio-LXLib (0.1.1)},
          year  = {2024},
          url   = {https://riskfolio-lib.readthedocs.io/en/latest/excel.html},
          }


License
=======

.. include:: ../../LICENSE-XL.txt